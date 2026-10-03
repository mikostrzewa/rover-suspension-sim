/*
 * Rover Suspension Sizer — physics module
 * =======================================
 * All math used by rover-suspension-sizer.html lives in this file. The page
 * only draws and collects inputs; every number it shows comes from a call into
 * RoverPhysics below. Nothing here touches the DOM, so the file also runs in
 * Node:   const RP = require('./physics.js');
 * Run the checks with:   node physics-tests.js   (or open verify.html)
 *
 * Sections
 *   1. Statics       2D level-body statics of one side of the rover (mm, N, deg)
 *   2. Terrain       ground profile: spline + obstacles, sampled polyline (mm / m)
 *   3. Vehicle       converts the design into SI parameters for the dynamics
 *   4. Dynamics      2D multibody equations of motion + RK4 integrator (SI)
 *   5. Scenarios     initial states for the terrain run and the drop test
 *   6. Tuning        spring rate / preload / damping search that minimises CG acceleration
 *
 * ---------------------------------------------------------------------------
 * Geometry and sign conventions (shared by statics and dynamics)
 * ---------------------------------------------------------------------------
 *   Origin at the body pivot O. x forward, y up. Pitch φ counter-clockwise
 *   (nose up for a rover driving towards +x).
 *   Two legs, both hinged at O:
 *       leg 0 = rear  = "Left β"  in the statics drawing, side sign s = −1
 *       leg 1 = front = "Right α" in the statics drawing, side sign s = +1
 *   θ = leg angle below the body's horizontal, measured outward.
 *   Body-frame points of a leg (L = leg length, a/e = mount 2 along / square to
 *   the leg, x1/y1 = mount 1 on the body, x outward):
 *       axle      b(θ)  = ( s·L·cosθ , −L·sinθ )
 *       mount 2   P2(θ) = ( s·(a·cosθ + e·sinθ) , −a·sinθ + e·cosθ )
 *       mount 1   P1    = ( s·x1 , y1 )
 *       spring    ℓ(θ)  = |P2 − P1|,  û = (P2 − P1)/ℓ   (points from mount 1 to 2)
 *   Coilover force F > 0 pushes the eyes apart (normal compression load).
 *   Wheel travel h > 0 is bump (wheel up relative to the body), h < 0 droop.
 *       sinθ(h) = sinθ0 − h/L
 */
(function (root, factory) {
  const RP = factory();
  if (typeof module === 'object' && module.exports) module.exports = RP;
  else root.RoverPhysics = RP;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const G = 9.81;                       // m/s²
  const D2R = Math.PI / 180, R2D = 180 / Math.PI;
  const clamp = (v, a, b) => Math.min(b, Math.max(a, v));

  /* =========================================================================
   * 1. STATICS  (units: mm, N, kg, degrees)
   * =========================================================================
   * Body held level. One side of the rover carries Fw = m·g·share at the CG.
   * Ground forces are vertical (free-rolling wheels, no traction).
   *
   * Whole side, ΣFy = 0 and ΣM_O = 0 (wheel contact directly below the axle):
   *     FN_front + FN_rear = Fw
   *     FN_front·X_front + FN_rear·X_rear = Fw·x_cg
   *   → FN_front = Fw·(x_cg − X_rear)/(X_front − X_rear)
   *
   * One leg, moments about O (pin reaction has no moment there):
   *     FN·x_w + Fs·(P2 × û) = 0          (× = 2D cross product)
   *   → Fs = FN · IR,   IR = −x_w / (P2 × û)          "leverage"
   *   Here x_w is the outward axle offset, so |P2 × û| is the perpendicular
   *   distance from O to the spring line (d⊥ in the sketch).
   * Pin reaction on the leg from ΣF = 0:   R = −(FN·ŷ + Fs·û)
   * Virtual work: Fs·dΔ = FN·dh, so the motion ratio MR = dΔ/dh = 1/IR.
   */

  /** Pose of one leg at wheel travel h. Leg in mm/deg: {L, th, x1, y1, a, e}.
   *  Returns outward-frame quantities (x positive towards this leg's wheel),
   *  or null when the leg cannot reach that travel. */
  function legPose(leg, h) {
    const s = Math.sin(leg.th * D2R) - h / leg.L;
    if (!(s >= -1 && s <= 1)) return null;
    const t = Math.asin(s), c = Math.cos(t), sn = Math.sin(t);
    const P2 = [leg.a * c + leg.e * sn, -leg.a * sn + leg.e * c];
    const W = [leg.L * c, -leg.L * sn];                    // axle
    const sv = [P2[0] - leg.x1, P2[1] - leg.y1];
    const l = Math.hypot(sv[0], sv[1]);
    if (l < 1e-6) return null;
    const u = [sv[0] / l, sv[1] / l];
    const arm = P2[0] * u[1] - P2[1] * u[0];               // P2 × û  (mm)
    return { t: t * R2D, P2, W, l, u, arm, IR: -W[0] / arm, MR: -arm / W[0] };
  }

  /** Level-body wheel loads (N). rover = {m [kg], share [%], cgx [mm]}. */
  function levelWheelLoads(rover, poseRear, poseFront) {
    const Wt = rover.m * G * rover.share / 100;
    const XF = poseFront.W[0], XR = -poseRear.W[0];
    const F = Wt * (rover.cgx - XR) / (XF - XR);
    return { Wt, L: Wt - F, R: F };
  }

  /** Force in a linear spring with preload, never negative (spring goes slack). */
  function springForceStatic(r, l) { return Math.max(0, r.k * (r.x0 + r.g0.l - l)); }
  function lengthAtTravel(leg, h) { const q = legPose(leg, h); return q ? q.l : NaN; }
  /** Bisection: wheel travel in [h0, h1] at which the spring length equals target. */
  function travelForLength(leg, target, h0, h1) {
    let a = h0, b = h1, fa = lengthAtTravel(leg, a) - target;
    for (let i = 0; i < 60; i++) {
      const m = (a + b) / 2, fm = lengthAtTravel(leg, m) - target;
      if (!Number.isFinite(fm)) return NaN;
      if ((fa < 0) === (fm < 0)) { a = m; fa = fm; } else b = m;
    }
    return (a + b) / 2;
  }

  /**
   * Complete statics + spring sizing for one side.
   * design = {m, share, cgx, mode, G, f, kman, legs:{L:{...}, R:{...}}}
   *   leg = {L, th, x1, y1, a, e, hmax, droop, ext, stroke} (mm, deg)
   * travel = {L: h_rear, R: h_front} in mm (pose used for the "current" values)
   *
   * Spring sizing (k in N/mm). Ride height is always the equilibrium anchor:
   *     F0 = IR(0)·FN0   and   preload x0 = F0/k
   * The second condition depends on mode:
   *   'droop'  spring free (F = 0) at full droop:  k = F0 / (ℓ_droop − ℓ0)
   *   'G'      wheel load G·FN0 at max bump:        k = (G·FN0·IR_bump − F0)/(ℓ0 − ℓ_bump)
   *   'freq'   ride frequency f:                    k = (2πf)²·(FN0/g)/MR0²   (/1000 → N/mm)
   *   'manual' k given
   * Wheel rate k_w = k·MR0², ride frequency f = √(k_w·1000 / m_c)/(2π), m_c = FN0/g.
   */
  function staticModel(design, travel) {
    const P = design.legs, out = { ok: true, msgs: [], r: {} };
    const G0 = {}, GN = {}, GB = {}, GD = {};
    for (const s of ['L', 'R']) {
      const p = P[s];
      G0[s] = legPose(p, 0); GN[s] = legPose(p, travel[s]); GB[s] = legPose(p, p.hmax); GD[s] = legPose(p, -p.droop);
      if (!G0[s] || !GN[s] || !GB[s] || !GD[s]) {
        out.ok = false;
        out.msgs.push(['bad', (s === 'L' ? 'Left' : 'Right') + ' leg cannot reach that travel. Reduce max bump or droop, or change the start angle or leg length.']);
      }
    }
    if (!out.ok) return out;
    const l0 = levelWheelLoads(design, G0.L, G0.R), ln = levelWheelLoads(design, GN.L, GN.R);
    out.Wt = l0.Wt;
    for (const s of ['L', 'R']) {
      const p = P[s], g0 = G0[s], gn = GN[s], gb = GB[s], gd = GD[s], nm = s === 'L' ? 'Left' : 'Right';
      const FN0 = l0[s], F0 = g0.IR * FN0, mc = FN0 / G, MR0 = g0.MR;
      const dB = g0.l - gb.l, dD = gd.l - g0.l;           // compression to bump, extension to droop (mm)
      let k;
      if (design.mode === 'G') k = (design.G * FN0 * gb.IR - F0) / dB;
      else if (design.mode === 'droop') k = F0 / dD;
      else if (design.mode === 'freq') k = Math.pow(2 * Math.PI * design.f, 2) * mc / 1000 / (MR0 * MR0);
      else k = design.kman;
      const r = { g0, gn, gb, gd, FN0, F0, mc, MR0, dB, dD, k };
      if (dB <= 0) out.msgs.push(['bad', nm + ' spring gets longer as the wheel moves up. Move mount 1 or mount 2 so the coilover compresses in bump.']);
      if (design.mode === 'droop' && dD <= 0) out.msgs.push(['bad', nm + ' spring does not extend in droop, so droop sizing cannot work.']);
      r.valid = Number.isFinite(k) && k > 0;
      if (!r.valid) {
        out.msgs.push(['bad', nm + ' spring rate comes out at zero or below with this method. ' + (design.mode === 'G' ? 'Raise the load at max bump.' : 'Check the geometry.')]);
        r.k = NaN;
      }
      r.x0 = F0 / k;                                       // preload compression at ride (mm)
      r.Fb = r.valid ? k * (r.x0 + dB) : NaN;              // spring force at max bump
      r.Gb = r.Fb / gb.IR / FN0;                           // wheel load at max bump ÷ static
      r.Fd = r.valid ? k * (r.x0 - dD) : NaN;              // spring force at full droop (<0 → slack)
      r.kw = k * MR0 * MR0;
      r.f = Math.sqrt(r.kw * 1000 / mc) / (2 * Math.PI);
      r.stroke = gd.l - gb.l;
      // current pose, 1 g statics
      const FN = ln[s], Fh = gn.IR * FN;
      r.FN = FN; r.Fh = Fh; r.pin = [-Fh * gn.u[0], -(FN + Fh * gn.u[1])]; r.comp = g0.l - gn.l;
      r.Fspr = r.valid ? springForceStatic(r, gn.l) : NaN; r.Fwh = r.Fspr / gn.IR;
      if (FN0 <= 0) out.msgs.push(['bad', 'The CG sits outside the wheelbase, so the ' + nm.toLowerCase() + ' wheel would lift off.']);
      if (gn.IR < 0) out.msgs.push(['bad', nm + ' spring would be pulled in tension at this pose. Its line of action passes the wrong side of the pivot.']);
      if (Math.abs(gn.arm) < 12) out.msgs.push(['warn', nm + ' spring line passes within 12 mm of the pivot, so the spring force is huge.']);
      if (r.valid && r.Fd < -0.5) out.msgs.push(['warn', nm + ' spring goes slack ' + Math.round(-r.Fd / k) + ' mm before full droop. Add preload or reduce droop.']);
      // coilover fit
      const c = { ext: p.ext, stroke: p.stroke, comp: p.ext - p.stroke, status: 'ok', text: 'Fits the whole travel' };
      if (p.stroke <= 0 || p.ext <= p.stroke) { c.status = 'bad'; c.text = 'Check the shock length and stroke'; }
      else if (g0.l < c.comp || g0.l > c.ext) { c.status = 'bad'; c.text = 'Ride length ' + g0.l.toFixed(1) + ' mm is outside the shock range'; }
      else {
        if (gb.l < c.comp) { c.status = 'bad'; c.hb = travelForLength(p, c.comp, 0, p.hmax); c.text = 'Bottoms out at +' + c.hb.toFixed(0) + ' mm bump (max bump +' + p.hmax + ')'; }
        if (gd.l > c.ext) {
          const hd = -travelForLength(p, c.ext, -p.droop, 0); c.hd = hd;
          const t = 'Tops out at −' + hd.toFixed(0) + ' mm droop (max droop −' + p.droop + ')';
          if (c.status === 'ok') { c.status = 'warn'; c.text = t; } else c.text += '; ' + t.toLowerCase();
        }
      }
      c.use = r.stroke / p.stroke; r.coil = c;
      if (c.status !== 'ok') out.msgs.push([c.status, nm + ' coilover: ' + c.text.charAt(0).toLowerCase() + c.text.slice(1) + '.']);
      out.r[s] = r;
    }
    return out;
  }

  /** Sweep one leg through its travel (other leg at ride). Returns rows of
   *  {h, l, comp, MR, FN, Fh (1 g hold force), Fs (sized spring), Fw (wheel load from Fs)}. */
  function travelSweep(design, M, s, N) {
    N = N || 60;
    const p = design.legs[s], o = s === 'L' ? 'R' : 'L', q0 = M.r[o].g0, r = M.r[s], rows = [];
    for (let i = 0; i <= N; i++) {
      const h = -p.droop + (p.hmax + p.droop) * i / N, gh = legPose(p, h);
      if (!gh) continue;
      const ld = s === 'L' ? levelWheelLoads(design, gh, q0) : levelWheelLoads(design, q0, gh);
      const FN = ld[s], Fh = gh.IR * FN, Fs = r.valid ? springForceStatic(r, gh.l) : NaN;
      rows.push({ h, l: gh.l, comp: r.g0.l - gh.l, MR: gh.MR, Fh, Fs, Fw: Fs / gh.IR, FN });
    }
    return rows;
  }

  /* =========================================================================
   * 2. TERRAIN
   * =========================================================================
   * Ground height (mm) = monotone cubic spline through control points
   *                      + sum of obstacle shape functions.
   * Sampled every 5 mm into a polyline; the dynamics uses that polyline in metres.
   * Obstacle shape f(d, o): d = x − o.x (mm), o = {h, w, n}.
   */
  const smoothstep = t => t <= 0 ? 0 : t >= 1 ? 1 : t * t * (3 - 2 * t);
  const OBSTACLES = {
    bump:     { name: 'Bump',      h: 50,  w: 400,  f: (d, o) => Math.abs(d) < o.w / 2 ? o.h * (1 + Math.cos(2 * Math.PI * d / o.w)) / 2 : 0, ext: o => o.w },
    rock:     { name: 'Rock',      h: 90,  w: 220,  f: (d, o) => { const q = Math.abs(2 * d / o.w); return q < 1 ? o.h * Math.pow(1 - q * q * q * q, 0.25) : 0; }, ext: o => o.w },
    stepup:   { name: 'Step up',   h: 80,  w: 10,   f: (d, o) => o.h * smoothstep((d + o.w / 2) / o.w), ext: o => o.w, end: true },
    stepdown: { name: 'Step down', h: 80,  w: 10,   f: (d, o) => -o.h * smoothstep((d + o.w / 2) / o.w), ext: o => o.w, end: true },
    ramp:     { name: 'Ramp',      h: 200, w: 1200, f: (d, o) => o.h * clamp((d + o.w / 2) / o.w, 0, 1), ext: o => o.w, end: true },
    ditch:    { name: 'Ditch',     h: 80,  w: 500,  f: (d, o) => Math.abs(d) < o.w / 2 ? -o.h * (1 + Math.cos(2 * Math.PI * d / o.w)) / 2 : 0, ext: o => o.w },
    wash:     { name: 'Washboard', h: 25,  w: 250, n: 6, f: (d, o) => { const L = o.n * o.w; return Math.abs(d) < L / 2 ? o.h * (1 - Math.cos(2 * Math.PI * (d + L / 2) / o.w)) / 2 : 0; }, ext: o => o.n * o.w }
  };

  function presetMap(name) {
    const ob = (type, x, h, w, n) => ({ type, x, h, w, n: n || 1 });
    switch (name) {
      case 'flat':  return { pts: [[0, 0], [10000, 0]], obs: [] };
      case 'rocks': return { pts: [[0, 0], [11000, 0]], obs: [ob('rock', 2400, 60, 180), ob('rock', 3500, 90, 240), ob('rock', 4300, 45, 150), ob('rock', 5200, 110, 300), ob('rock', 6500, 70, 200), ob('rock', 7300, 120, 280), ob('rock', 8700, 80, 220)] };
      case 'steps': return { pts: [[0, 0], [11000, 0]], obs: [ob('stepup', 2400, 50, 10), ob('stepup', 3800, 80, 10), ob('stepup', 5200, 100, 10), ob('stepdown', 7000, 100, 10), ob('stepdown', 8600, 130, 10)] };
      case 'wash':  return { pts: [[0, 0], [12000, 0]], obs: [ob('wash', 6000, 25, 300, 24)] };
      case 'slope': { const L = 3000, hh = Math.round(L * Math.tan(15 * D2R)); return { pts: [[0, 0], [2500, 0], [2500 + L, hh], [2500 + L + 3000, hh]], obs: [] }; }
      default:      return { pts: [[0, 0], [5000, 0], [6500, 250], [9600, 250], [11000, 0], [15500, 0]], obs: [ob('bump', 2300, 50, 400), ob('rock', 3800, 90, 220), ob('stepdown', 8400, 100, 10), ob('ditch', 11800, 80, 500), ob('wash', 13200, 25, 250, 6)] };
    }
  }

  /** Fritsch–Carlson monotone cubic Hermite spline through sorted [x, y] points.
   *  Flat beyond the end points. Never overshoots between control points. */
  function monotoneSpline(pts) {
    const n = pts.length, xs = pts.map(p => p[0]), ys = pts.map(p => p[1]);
    if (n < 2) return () => ys[0] || 0;
    const d = [], m = new Array(n);
    for (let k = 0; k < n - 1; k++) d.push((ys[k + 1] - ys[k]) / Math.max(1e-6, xs[k + 1] - xs[k]));
    m[0] = d[0]; m[n - 1] = d[n - 2];
    for (let k = 1; k < n - 1; k++) m[k] = d[k - 1] * d[k] <= 0 ? 0 : (d[k - 1] + d[k]) / 2;
    for (let k = 0; k < n - 1; k++) {
      if (d[k] === 0) { m[k] = 0; m[k + 1] = 0; continue; }
      const a = m[k] / d[k], b = m[k + 1] / d[k], h = a * a + b * b;
      if (h > 9) { const t = 3 / Math.sqrt(h); m[k] = t * a * d[k]; m[k + 1] = t * b * d[k]; }
    }
    return x => {
      if (x <= xs[0]) return ys[0];
      if (x >= xs[n - 1]) return ys[n - 1];
      let lo = 0, hi = n - 1;
      while (hi - lo > 1) { const md = (lo + hi) >> 1; if (xs[md] <= x) lo = md; else hi = md; }
      const h = xs[hi] - xs[lo], t = (x - xs[lo]) / h, t2 = t * t, t3 = t2 * t;
      return (2 * t3 - 3 * t2 + 1) * ys[lo] + (t3 - 2 * t2 + t) * h * m[lo] + (-2 * t3 + 3 * t2) * ys[hi] + (t3 - t2) * h * m[hi];
    };
  }

  /** Build the sampled terrain for a map {pts:[[x,y] mm], obs:[...]}.
   *  Returns mm-space helpers for drawing and `ter` (metres) for the dynamics. */
  function buildTerrain(map, dx) {
    dx = dx || 5;
    const pts = map.pts, spl = monotoneSpline(pts), len = pts[pts.length - 1][0], n = Math.floor(len / dx) + 1;
    const Y = new Float64Array(n), Ym = new Float64Array(n);
    let lo = Infinity, hi = -Infinity;
    const ground = x => { let y = spl(x); for (const o of map.obs) y += OBSTACLES[o.type].f(x - o.x, o); return y; };
    for (let i = 0; i < n; i++) { const y = ground(i * dx); Y[i] = y; Ym[i] = y / 1000; if (y < lo) lo = y; if (y > hi) hi = y; }
    const hmm = x => { if (x <= 0) return Y[0]; if (x >= (n - 1) * dx) return Y[n - 1]; const i = Math.floor(x / dx), f = x / dx - i; return Y[i] * (1 - f) + Y[i + 1] * f; };
    const groundWithout = (x, skip) => { let y = spl(x); for (const o of map.obs) if (o !== skip) y += OBSTACLES[o.type].f(x - o.x, o); return y; };
    return { spl, len, dx, n, Y, lo, hi, hmm, ground, groundWithout, ter: contactTerrain(Ym, dx / 1000) };
  }

  /**
   * Contact polyline in metres: Y[i] is the height at x = i·dx.
   * contacts(px, py, r, out) finds where a wheel (centre px,py, radius r)
   * overlaps the ground. For each polyline segment near the wheel it takes the
   * closest point; a contact patch is a local minimum of that distance along the
   * polyline with distance < r (so a flat floor gives one contact, a step edge
   * plus floor gives two). Writes up to 3 contacts as [x, y, distance] triples.
   */
  function contactTerrain(Y, dx) {
    const n = Y.length;
    const T = {
      dx, n, Y, len: (n - 1) * dx,
      h(x) { if (x <= 0) return Y[0]; if (x >= T.len) return Y[n - 1]; const i = Math.floor(x / dx), f = x / dx - i; return Y[i] * (1 - f) + Y[i + 1] * f; },
      cx: new Float64Array(160), cy: new Float64Array(160), cd: new Float64Array(160),
      contacts(px, py, r, out) {
        const i0 = Math.max(0, Math.floor((px - r) / dx) - 1), i1 = Math.min(n - 2, Math.ceil((px + r) / dx) + 1);
        let m = 0;
        for (let i = i0; i <= i1 && m < 160; i++) {
          const x1 = i * dx, y1 = Y[i], ex = dx, ey = Y[i + 1] - y1;
          let t = ((px - x1) * ex + (py - y1) * ey) / (ex * ex + ey * ey);
          t = t < 0 ? 0 : t > 1 ? 1 : t;
          const qx = x1 + t * ex, qy = y1 + t * ey;
          T.cx[m] = qx; T.cy[m] = qy; T.cd[m] = Math.hypot(px - qx, py - qy); m++;
        }
        let k = 0;
        for (let j = 0; j < m && k < 3; j++) {
          const dj = T.cd[j]; if (dj >= r) continue;
          const L = j > 0 ? T.cd[j - 1] : Infinity, R = j < m - 1 ? T.cd[j + 1] : Infinity;
          if (dj <= L && dj < R) { out[k * 3] = T.cx[j]; out[k * 3 + 1] = T.cy[j]; out[k * 3 + 2] = dj; k++; }
        }
        return k;
      }
    };
    return T;
  }
  let FLAT = null;
  /** 30 m of flat ground at y = 0 (drop test). */
  const flatTerrain = () => FLAT || (FLAT = contactTerrain(new Float64Array(3001), 0.01));

  /* =========================================================================
   * 3. VEHICLE  (design in mm/kg → SI vehicle for the dynamics)
   * =========================================================================
   * Mass split of one side (m_side = m·share):
   *   body (sprung) m_b at the CG, pitch inertia I_b = m_b·k_y²
   *   per leg: wheel-end mass m_w at the axle (wheel + motor), leg as a uniform
   *   rod of mass m_l from O to the axle (I_rod = m_l·L²/12 about its middle)
   *   m_b = m_side − 2(m_w + m_l)
   * Spring rate k comes from the statics sizing. Preload is re-solved for this
   * mass split so the rover rests exactly at θ0 on flat ground:
   *   static wheel loads from the total centroid
   *   x̄ = (m_b·x_cg + (m_w + m_l/2)(X_rear + X_front)) / m_side
   *   leg moments: F_spring = IR·(FN − (m_w + m_l/2)·g),  x0 = F_spring / k
   * Tyre: 'rigid' → stiff penalty contact k_t = 2·10⁶ N/m with damping for
   * ζ = 0.9 on the wheel-end mass (≈0.06 mm deflection under 120 N);
   * 'compliant' → the k_t and c_t you enter.
   */
  const RIGID_KT = 2e6, RIGID_ZETA = 0.9;
  function buildVehicle(design, M) {
    if (!M || !M.ok) return { err: 'Fix the geometry problem on the Statics tab first.' };
    const D = design.dyn, side = design.m * design.share / 100, mw = D.mw, ml = D.mleg, mb = side - 2 * (mw + ml);
    if (!(mb > 0.5)) return { err: 'Wheel-end and leg masses use up nearly all the weight on this side. Lower them or raise the rover mass.' };
    const legs = [], Xw = [], ids = ['L', 'R'];
    for (let i = 0; i < 2; i++) {
      const s = ids[i], p = design.legs[s], r = M.r[s];
      if (!r.valid) return { err: 'The spring rate on the Statics tab is not valid, so there is no spring to simulate.' };
      legs.push({ s: i ? 1 : -1, L: p.L / 1000, a: p.a / 1000, e: p.e / 1000, x1: p.x1 / 1000, y1: p.y1 / 1000,
        k: r.k * 1000, x0: 0, l0: r.g0.l / 1000, lmin: (p.ext - p.stroke) / 1000, lmax: p.ext / 1000,
        th0: p.th * D2R, IR: r.g0.IR, MR0: r.MR0 });
      Xw.push((i ? 1 : -1) * p.L / 1000 * Math.cos(p.th * D2R));
    }
    const xt = (mb * design.cgx / 1000 + (mw + ml / 2) * (Xw[0] + Xw[1])) / side;
    const FNF = side * G * (xt - Xw[0]) / (Xw[1] - Xw[0]), FN = [side * G - FNF, FNF];
    if (FN[0] <= 0 || FN[1] <= 0) return { err: 'The CG sits outside the wheelbase, so one wheel would lift off at rest.' };
    for (let i = 0; i < 2; i++) { const Fs = legs[i].IR * (FN[i] - (mw + ml / 2) * G); legs[i].x0 = Fs / legs[i].k; legs[i].F0 = Fs; }
    const rigid = D.tyre !== 'compliant';
    const kt = rigid ? RIGID_KT : D.kt * 1000;
    const ct = rigid ? 2 * RIGID_ZETA * Math.sqrt(RIGID_KT * mw) : D.ct;
    return {
      mb, Ib: mb * Math.pow(D.ky / 1000, 2), mw, ml, r: design.rw / 1000, J: Math.max(1e-4, D.J),
      kt, ct, rigid, mu: D.mu, crr: D.crr, cb: D.cb, cr: D.cr, kbs: D.kbs * 1000, cbs: 2500,
      Ts: D.Ts, w0: Math.max(0.1, D.rpm) * 2 * Math.PI / 60, drive: [D.dRear, D.dFront], coast: D.coast,
      cg: [design.cgx / 1000, design.cgy / 1000], legs, g: G, vEps: 0.05, FN, side
    };
  }

  /* =========================================================================
   * 4. DYNAMICS  (SI units)
   * =========================================================================
   * State vector st (14):
   *   0 x   1 y   2 φ   3 θ_rear   4 θ_front          generalized coordinates q
   *   5..9  q̇
   *   10 ω_rear  11 ω_front     absolute wheel spin, + = clockwise = forward roll
   *   12 ψ_rear  13 ψ_front     wheel angles (drawing only)
   *
   * Kinetic energy, R(φ) = rotation, β_i = φ − s_i·θ_i = absolute leg angle:
   *   T = ½m_b|v_cg|² + ½I_b·φ̇² + Σ_i [ ½m_w|v_axle,i|² + ½m_l|v_mid,i|²
   *                                    + ½I_rod·β̇_i² + ½J·ω_i² ]
   * Every mass point p = O + R(φ)·b(q) has velocity v = J_p·q̇ with columns
   *   ∂p/∂x = (1,0), ∂p/∂y = (0,1), ∂p/∂φ = R·S·b, ∂p/∂θ_i = R·b'(θ_i)
   * (S = 90° rotation, so R·S·b = perp(R·b)). Lagrange's equations become
   *   M(q)·q̈ = Q − Σ_k m_k·J_kᵀ·a_vp,k
   *   M = Σ m_k·J_kᵀJ_k + I_b·e_φe_φᵀ + Σ I_rod·∇β_i∇β_iᵀ
   *   a_vp = velocity-product (centripetal/Coriolis) part of the acceleration:
   *     axle:  R·( −b·φ̇² + 2·S·b'·φ̇·θ̇ − b·θ̇² )      (uses b'' = −b)
   *     CG:    −R·c·φ̇²
   * Generalized forces Q:
   *   gravity              J_kᵀ·(0, −m_k·g)
   *   coilover             F·∂ℓ/∂θ_i   in the θ_i row (spring + damper + stops)
   *   tyre contact         J_axleᵀ·F_c        (F_c acts on the wheel; the axle carries it)
   *   motor reaction       +T_i·∇β_i          (motor on the leg; wheel gets −T)
   * Wheel spin:  J·ω̇_i = T_i − τ_c,i − M_rr,i
   *   τ_c = Σ (p_c − p_axle) × F_c  (counter-clockwise torque of contact forces)
   *   M_rr = C_rr·N·r·tanh(ω/0.5)    rolling resistance moment
   * Tyre contact (penalty): δ = r − |p_axle − p_c|,  n̂ from contact to axle,
   *   N = max(0, k_t·δ + c_t·δ̇),  slip s = (v_rim at contact)·t̂,
   *   F_t = −μ·N·tanh(s / v_ε),  F_c = N·n̂ + F_t·t̂,  t̂ = (n_y, −n_x)
   * Motor (DC curve, per wheel):  T = clamp(T_stall·(u − ω_rel/ω_0), ±T_stall),
   *   ω_rel = ω + β̇ (wheel speed relative to the leg the motor is bolted to).
   * Coilover:  F = F_spring + F_damper + F_stop
   *   F_spring = max(0, k·(x0 + ℓ0 − ℓ))
   *   F_damper = −c_bump·ℓ̇ (ℓ̇ < 0)  or  −c_rebound·ℓ̇ (ℓ̇ > 0)
   *   F_stop   = k_bs·(ℓ_min − ℓ) + c·(compression) below ℓ_min, mirrored above ℓ_max
   */

  /** Body-frame leg kinematics at angle th (rad), leg in metres (from buildVehicle). */
  function legKinematics(lg, th) {
    const c = Math.cos(th), sn = Math.sin(th), s = lg.s;
    const bx = s * lg.L * c, by = -lg.L * sn;                    // axle b(θ)
    const dbx = -s * lg.L * sn, dby = -lg.L * c;                 // b'(θ)
    const p2x = s * (lg.a * c + lg.e * sn), p2y = -lg.a * sn + lg.e * c;     // mount 2
    const d2x = s * (-lg.a * sn + lg.e * c), d2y = -lg.a * c - lg.e * sn;   // dP2/dθ
    const dx = p2x - s * lg.x1, dy = p2y - lg.y1, l = Math.hypot(dx, dy);
    const dl = (dx * d2x + dy * d2y) / l;                         // ∂ℓ/∂θ
    return { bx, by, dbx, dby, p2x, p2y, l, dl };
  }

  /** Coilover force along the shock (N, + pushes eyes apart). Optional out gets the parts. */
  function coiloverForce(V, lg, l, ldot, out) {
    const Fs = Math.max(0, lg.k * (lg.x0 + lg.l0 - l));
    const Fd = ldot < 0 ? -V.cb * ldot : -V.cr * ldot;
    let Fb = 0;
    if (l < lg.lmin) Fb = V.kbs * (lg.lmin - l) + (ldot < 0 ? -V.cbs * ldot : 0);
    else if (l > lg.lmax) Fb = -V.kbs * (l - lg.lmax) - (ldot > 0 ? V.cbs * ldot : 0);
    if (out) { out.Fs = Fs; out.Fd = Fd; out.Fb = Fb; }
    return Fs + Fd + Fb;
  }

  /** In-place Gaussian elimination with partial pivoting: solves A·x = b (A is n×n row-major). */
  function solveLinear(A, b, n) {
    for (let c = 0; c < n; c++) {
      let p = c, mx = Math.abs(A[c * n + c]);
      for (let r = c + 1; r < n; r++) { const v = Math.abs(A[r * n + c]); if (v > mx) { mx = v; p = r; } }
      if (p !== c) { for (let k = 0; k < n; k++) { const t = A[c * n + k]; A[c * n + k] = A[p * n + k]; A[p * n + k] = t; } const t = b[c]; b[c] = b[p]; b[p] = t; }
      const d = A[c * n + c];
      for (let r = c + 1; r < n; r++) { const f = A[r * n + c] / d; if (f === 0) continue; for (let k = c; k < n; k++) A[r * n + k] -= f * A[c * n + k]; b[r] -= f * b[c]; }
    }
    for (let r = n - 1; r >= 0; r--) { let s = b[r]; for (let k = r + 1; k < n; k++) s -= A[r * n + k] * b[k]; b[r] = s / A[r * n + r]; }
    return b;
  }

  const cbuf = new Float64Array(9);
  /** Allocate the work arrays a vehicle needs (called once by Simulation). */
  function prepareVehicle(V) {
    V._M = new Float64Array(25); V._Q = new Float64Array(5); V._Jx = new Float64Array(5); V._Jy = new Float64Array(5); V.Nmax = 0;
    V.legs.forEach(lg => { lg.Irod = V.ml * lg.L * lg.L / 12; });
    return V;
  }

  /**
   * Assemble M (5×5, into V._M) and the right-hand side Q − Σ m·Jᵀ·a_vp (into V._Q),
   * and write the wheel-spin accelerations into d[10..13].
   * u = throttle in [−1, 1]. inf (optional) receives forces for recording.
   */
  function assemble(V, T, st, u, d, inf) {
    const g = V.g, ph = st[2], vx = st[5], vy = st[6], vp = st[7], c = Math.cos(ph), s = Math.sin(ph);
    const M = V._M, Q = V._Q, Jx = V._Jx, Jy = V._Jy;
    M.fill(0); Q.fill(0);
    // point mass m with Jacobian rows Jx, Jy and velocity-product acceleration (ax, ay)
    function addPoint(m, ax, ay) {
      for (let i = 0; i < 5; i++) {
        const xi = Jx[i], yi = Jy[i];
        if (xi === 0 && yi === 0) continue;
        for (let j = 0; j < 5; j++) M[i * 5 + j] += m * (xi * Jx[j] + yi * Jy[j]);
        Q[i] += -m * (xi * ax + yi * (ay + g));         // −m·Jᵀ·a_vp  plus gravity Jᵀ·(0, −m·g)
      }
    }
    // body: CG at c (body frame), rotated Rc
    const rcx = c * V.cg[0] - s * V.cg[1], rcy = s * V.cg[0] + c * V.cg[1];
    Jx.fill(0); Jy.fill(0); Jx[0] = 1; Jy[1] = 1; Jx[2] = -rcy; Jy[2] = rcx;
    const acx = -rcx * vp * vp, acy = -rcy * vp * vp;
    addPoint(V.mb, acx, acy);
    M[2 * 5 + 2] += V.Ib;
    let Nmax = 0;
    for (let i = 0; i < 2; i++) {
      const lg = V.legs[i], k = 3 + i, th = st[3 + i], vt = st[8 + i], w = st[10 + i];
      const K = legKinematics(lg, th);
      const rbx = c * K.bx - s * K.by, rby = s * K.bx + c * K.by;       // R·b
      const rdx = c * K.dbx - s * K.dby, rdy = s * K.dbx + c * K.dby;   // R·b'
      // a_vp of the axle: R(−b φ̇² + 2 S b' φ̇ θ̇ − b θ̇²);  S·v = (−v_y, v_x)
      const avx = -rbx * vp * vp + 2 * (-rdy) * vp * vt - rbx * vt * vt;
      const avy = -rby * vp * vp + 2 * (rdx) * vp * vt - rby * vt * vt;
      Jx.fill(0); Jy.fill(0); Jx[0] = 1; Jy[1] = 1; Jx[2] = -rby; Jy[2] = rbx; Jx[k] = rdx; Jy[k] = rdy;
      addPoint(V.mw, avx, avy);                                         // wheel-end mass at the axle
      Jx[2] = -rby / 2; Jy[2] = rbx / 2; Jx[k] = rdx / 2; Jy[k] = rdy / 2;
      addPoint(V.ml, avx / 2, avy / 2);                                 // leg rod centre (b/2)
      const gs = -lg.s;                                                 // ∂β/∂θ_i
      M[2 * 5 + 2] += lg.Irod; M[2 * 5 + k] += lg.Irod * gs; M[k * 5 + 2] += lg.Irod * gs; M[k * 5 + k] += lg.Irod;
      // coilover: generalized force F·∂ℓ/∂θ
      const ld = K.dl * vt, co = inf ? inf.co[i] : null;
      const F = coiloverForce(V, lg, K.l, ld, co);
      Q[k] += F * K.dl;
      // tyre contact
      const pax = st[0] + rbx, pay = st[1] + rby, vax = vx - rby * vp + rdx * vt, vay = vy + rbx * vp + rdy * vt;
      const nc = T.contacts(pax, pay, V.r, cbuf);
      let Nt = 0, Ftt = 0, tau = 0, fxs = 0, fys = 0, slip0 = 0;
      for (let j = 0; j < nc; j++) {
        const qx = cbuf[j * 3], qy = cbuf[j * 3 + 1], dd = Math.max(1e-9, cbuf[j * 3 + 2]);
        const nx = (pax - qx) / dd, ny = (pay - qy) / dd;               // normal, contact → axle
        const del = V.r - dd, deld = -(vax * nx + vay * ny);            // penetration and its rate
        const N = Math.max(0, V.kt * del + V.ct * deld);
        if (N <= 0) continue;
        const tx = ny, ty = -nx, rx = qx - pax, ry = qy - pay;
        // rim velocity at contact: v_axle + ω_ccw × r_c,  ω_ccw = −ω
        const vcx = vax + (-w) * (-ry), vcy = vay + (-w) * (rx);
        const sl = vcx * tx + vcy * ty;
        const Ft = -V.mu * N * Math.tanh(sl / V.vEps);
        const Fx = N * nx + Ft * tx, Fy = N * ny + Ft * ty;
        Q[0] += Fx; Q[1] += Fy; Q[2] += (-rby) * Fx + rbx * Fy; Q[k] += rdx * Fx + rdy * Fy;   // J_axleᵀ·F
        tau += rx * Fy - ry * Fx;                                       // ccw torque about the axle
        Nt += N; Ftt += Ft; fxs += Fx; fys += Fy;
        if (j === 0) slip0 = sl;
        if (inf && j === 0) { inf.cp[i * 2] = qx; inf.cp[i * 2 + 1] = qy; }
      }
      if (Nt > Nmax) Nmax = Nt;
      // motor on the leg
      const bd = vp - lg.s * vt, wr = w + bd;
      let Tm = 0;
      if (V.drive[i] && !(V.coast && u === 0)) Tm = clamp(V.Ts * (u - wr / V.w0), -V.Ts, V.Ts);
      Q[2] += Tm; Q[k] += -lg.s * Tm;                                   // reaction T·∇β on the leg
      const Mrr = V.crr * Nt * V.r * Math.tanh(w / 0.5);
      d[10 + i] = (Tm - tau - Mrr) / V.J;
      d[12 + i] = w;
      if (inf) { inf.N[i] = Nt; inf.Ft[i] = Ftt; inf.F[i] = F; inf.l[i] = K.l; inf.T[i] = Tm; inf.slip[i] = slip0; inf.nc[i] = nc; inf.fx[i] = fxs; inf.fy[i] = fys; inf.pa[i * 2] = pax; inf.pa[i * 2 + 1] = pay; }
    }
    V.Nmax = Nmax;
    return { rcx, rcy, acx, acy };
  }

  /** State derivative d = f(st, u). */
  function derivatives(V, T, st, u, d, inf) {
    const a = assemble(V, T, st, u, d, inf);
    const qdd = solveLinear(V._M, V._Q, 5);
    for (let i = 0; i < 5; i++) { d[i] = st[5 + i]; d[5 + i] = qdd[i]; }
    if (inf) {   // accelerometer at the CG (load factor, g units)
      const ax = qdd[0] - a.rcy * qdd[2] + a.acx, ay = qdd[1] + a.rcx * qdd[2] + a.acy;
      inf.nz = (ay + V.g) / V.g; inf.nx = ax / V.g;
    }
    return d;
  }

  /** Copies of M and of the generalized-force vector (Q − Σ m·Jᵀ·a_vp) at a state. For verification. */
  function massMatrixAndForces(V, T, st, u) {
    if (!V._M) prepareVehicle(V);
    const d = new Float64Array(14);
    assemble(V, T, st, u || 0, d);
    return { M: Array.from(V._M), Q: Array.from(V._Q), wheelAcc: [d[10], d[11]] };
  }

  /** Total mechanical energy: kinetic + gravity + spring + stops + tyre. For verification. */
  function totalEnergy(V, T, st) {
    if (!V._M) prepareVehicle(V);
    const c = Math.cos(st[2]), s = Math.sin(st[2]);
    const rcx = c * V.cg[0] - s * V.cg[1], rcy = s * V.cg[0] + c * V.cg[1];
    let KE = 0.5 * V.mb * ((st[5] - rcy * st[7]) ** 2 + (st[6] + rcx * st[7]) ** 2) + 0.5 * V.Ib * st[7] ** 2;
    let PE = V.mb * V.g * (st[1] + rcy);
    for (let i = 0; i < 2; i++) {
      const lg = V.legs[i], K = legKinematics(lg, st[3 + i]), vt = st[8 + i];
      const rbx = c * K.bx - s * K.by, rby = s * K.bx + c * K.by, rdx = c * K.dbx - s * K.dby, rdy = s * K.dbx + c * K.dby;
      const vax = st[5] - rby * st[7] + rdx * vt, vay = st[6] + rbx * st[7] + rdy * vt;
      const vmx = st[5] + (-rby * st[7] + rdx * vt) / 2, vmy = st[6] + (rbx * st[7] + rdy * vt) / 2;
      KE += 0.5 * V.mw * (vax * vax + vay * vay) + 0.5 * V.ml * (vmx * vmx + vmy * vmy) + 0.5 * lg.Irod * (st[7] - lg.s * vt) ** 2 + 0.5 * V.J * st[10 + i] ** 2;
      PE += V.mw * V.g * (st[1] + rby) + V.ml * V.g * (st[1] + rby / 2);
      const x = lg.x0 + lg.l0 - K.l; if (x > 0) PE += 0.5 * lg.k * x * x;
      if (K.l < lg.lmin) PE += 0.5 * V.kbs * (lg.lmin - K.l) ** 2;
      if (K.l > lg.lmax) PE += 0.5 * V.kbs * (K.l - lg.lmax) ** 2;
      const nc = T.contacts(st[0] + rbx, st[1] + rby, V.r, cbuf);
      for (let j = 0; j < nc; j++) { const del = V.r - cbuf[j * 3 + 2]; PE += 0.5 * V.kt * del * del; }
    }
    return { KE, PE, E: KE + PE };
  }

  function newInfo() {
    return { co: [{ Fs: 0, Fd: 0, Fb: 0 }, { Fs: 0, Fd: 0, Fb: 0 }], N: [0, 0], Ft: [0, 0], F: [0, 0], l: [0, 0], T: [0, 0], slip: [0, 0], nc: [0, 0], fx: [0, 0], fy: [0, 0], cp: [0, 0, 0, 0], pa: [0, 0, 0, 0], nz: 1, nx: 0 };
  }
  /** Recorded channels (index 0 = rear, 1 = front). Units SI.
   *  Most channels are instantaneous values at the sample time. The *a channels
   *  (N0a, N1a, F0a, F1a, nza, nxa) are time averages over the interval since the
   *  previous sample, integrated at every RK4 step. They keep the impulse of very
   *  short spikes (a rigid wheel hitting a rock lasts ~1 ms) that point samples
   *  would miss or exaggerate, like a data logger with an anti-alias filter. */
  const RECORD_KEYS = ['t', 'x', 'y', 'ph', 'th0', 'th1', 'ps0', 'ps1', 'w0', 'w1', 'u', 'v', 'nz', 'nx', 'N0', 'N1', 'Ft0', 'Ft1', 'l0', 'l1', 'F0', 'F1',
    'Fs0', 'Fs1', 'Fd0', 'Fd1', 'Fb0', 'Fb1', 'T0', 'T1', 'sl0', 'sl1', 'cpx0', 'cpy0', 'cpx1', 'cpy1', 'fx0', 'fy0', 'fx1', 'fy1', 'nc0', 'nc1', 'P',
    'N0a', 'N1a', 'F0a', 'F1a', 'nza', 'nxa'];

  /**
   * Fixed-throttle RK4 simulation with an adaptive step.
   * Step size keeps every stiff term inside RK4's stability limit:
   *   friction   λ_f = μ·N_max/v_ε·(r²/J + 1/m_w)
   *   tyre       ω_t = √(k_t/m_w),  λ_c = c_t/m_w
   *   h = clamp(min(0.85/λ_f, 1/ω_t, 1/λ_c), 15 µs, dtMax)
   * Records every recDt seconds into this.rec (arrays keyed by RECORD_KEYS).
   */
  function Simulation(V, T, st0, opt) {
    this.V = prepareVehicle(V); this.T = T; this.st = Float64Array.from(st0); this.t = 0;
    this.k1 = new Float64Array(14); this.k2 = new Float64Array(14); this.k3 = new Float64Array(14); this.k4 = new Float64Array(14); this.tmp = new Float64Array(14);
    this.inf = newInfo(); this.rec = {}; RECORD_KEYS.forEach(k => this.rec[k] = []);
    this.recDt = (opt && opt.recDt) || 0.005; this.nextRec = 0; this.dtMax = (opt && opt.dtMax) || 4e-4; this.status = '';
    this.hTyre = Math.min(1 / Math.sqrt(V.kt / V.mw), V.ct > 0 ? 1 / (V.ct / V.mw) : Infinity);
    this.infK = [newInfo(), newInfo(), newInfo(), newInfo()]; this.acc = { t: 0, N: [0, 0], F: [0, 0], nz: 0, nx: 0 };
  }
  Simulation.prototype.dt = function () {
    const V = this.V, lam = V.mu * Math.max(V.Nmax, 30) / V.vEps * (V.r * V.r / V.J + 1 / V.mw) + 1e-9;
    return clamp(Math.min(0.85 / lam, this.hTyre), 1.5e-5, this.dtMax);
  };
  Simulation.prototype.step = function (h, u) {
    const V = this.V, T = this.T, y = this.st, n = 14, k1 = this.k1, k2 = this.k2, k3 = this.k3, k4 = this.k4, tp = this.tmp;
    const I = this.infK, A = this.acc;
    derivatives(V, T, y, u, k1, I[0]);
    for (let i = 0; i < n; i++) tp[i] = y[i] + h / 2 * k1[i]; derivatives(V, T, tp, u, k2, I[1]);
    for (let i = 0; i < n; i++) tp[i] = y[i] + h / 2 * k2[i]; derivatives(V, T, tp, u, k3, I[2]);
    for (let i = 0; i < n; i++) tp[i] = y[i] + h * k3[i]; derivatives(V, T, tp, u, k4, I[3]);
    for (let i = 0; i < n; i++) y[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    // interval averages with the same RK4 weights (consistent quadrature of the forces)
    const w6 = h / 6;
    for (let j = 0; j < 4; j++) { const wj = (j === 1 || j === 2 ? 2 : 1) * w6, Ij = I[j]; A.N[0] += Ij.N[0] * wj; A.N[1] += Ij.N[1] * wj; A.F[0] += Ij.F[0] * wj; A.F[1] += Ij.F[1] * wj; A.nz += Ij.nz * wj; A.nx += Ij.nx * wj; }
    A.t += h;
    this.t += h;
  };
  Simulation.prototype.record = function (u) {
    const I = this.inf, y = this.st, R = this.rec;
    derivatives(this.V, this.T, y, u, this.k4, I);
    const P = I.T[0] * (y[10] + y[7] - this.V.legs[0].s * y[8]) + I.T[1] * (y[11] + y[7] - this.V.legs[1].s * y[9]);   // motor power
    const v = [this.t, y[0], y[1], y[2], y[3], y[4], y[12], y[13], y[10], y[11], u, y[5], I.nz, I.nx, I.N[0], I.N[1], I.Ft[0], I.Ft[1], I.l[0], I.l[1], I.F[0], I.F[1],
      I.co[0].Fs, I.co[1].Fs, I.co[0].Fd, I.co[1].Fd, I.co[0].Fb, I.co[1].Fb, I.T[0], I.T[1], I.slip[0], I.slip[1],
      I.nc[0] ? I.cp[0] : NaN, I.nc[0] ? I.cp[1] : NaN, I.nc[1] ? I.cp[2] : NaN, I.nc[1] ? I.cp[3] : NaN, I.fx[0], I.fy[0], I.fx[1], I.fy[1], I.nc[0], I.nc[1], P];
    const A = this.acc, at = A.t > 0 ? A.t : 0;
    v.push(at ? A.N[0] / at : I.N[0], at ? A.N[1] / at : I.N[1], at ? A.F[0] / at : I.F[0], at ? A.F[1] / at : I.F[1], at ? A.nz / at : I.nz, at ? A.nx / at : I.nx);
    A.t = 0; A.N[0] = A.N[1] = A.F[0] = A.F[1] = A.nz = A.nx = 0;
    for (let i = 0; i < RECORD_KEYS.length; i++) R[RECORD_KEYS[i]].push(v[i]);
  };
  /** Integrate to tEnd with throttle u. Returns false if the state stops being finite. */
  Simulation.prototype.advance = function (tEnd, u) {
    while (this.t < tEnd - 1e-12) {
      if (this.t >= this.nextRec - 1e-12) { this.record(u); this.nextRec += this.recDt; }
      const h = Math.min(this.dt(), tEnd - this.t, Math.max(1e-6, this.nextRec - this.t));
      this.step(h, u);
      for (let i = 0; i < 14; i++) if (!Number.isFinite(this.st[i])) { this.status = 'diverged'; return false; }
    }
    return true;
  };

  /* =========================================================================
   * 5. SCENARIOS
   * ========================================================================= */
  /** World x of both axles. */
  function wheelXs(V, st) {
    const c = Math.cos(st[2]), s = Math.sin(st[2]);
    return [0, 1].map(i => { const K = legKinematics(V.legs[i], st[3 + i]); return st[0] + c * K.bx - s * K.by; });
  }
  /** Terrain run start: legs at θ0, body level, tyres pressed in by their static load. */
  function initialRunState(V, T, x0, v0) {
    const st = new Float64Array(14); st[0] = x0; let y = -1e9;
    for (let i = 0; i < 2; i++) {
      const lg = V.legs[i], K = legKinematics(lg, lg.th0); st[3 + i] = lg.th0;
      y = Math.max(y, T.h(x0 + K.bx) + V.r - V.FN[i] / V.kt - K.by);
    }
    st[1] = y; st[5] = v0; st[10] = st[11] = v0 / V.r;
    return st;
  }
  /** Leg angle hanging in the air (body held at pitch ph): coilover moment balances
   *  the leg's own weight,  F·∂ℓ/∂θ − (m_w + m_l/2)·g·(R·b')_y = 0. */
  function airborneLegAngle(V, lg, ph) {
    const c = Math.cos(ph), s = Math.sin(ph), mq = V.mw + V.ml / 2;
    const q = th => { const K = legKinematics(lg, th), F = coiloverForce(V, lg, K.l, 0); return F * K.dl - mq * V.g * (s * K.dbx + c * K.dby); };
    const a = lg.th0;
    if (q(a) <= 0) return a;
    for (let th = a + 0.005; th < 1.45; th += 0.005) {
      if (q(th) <= 0) { let lo = th - 0.005, hi = th; for (let k = 0; k < 40; k++) { const m = (lo + hi) / 2; if (q(m) > 0) lo = m; else hi = m; } return (lo + hi) / 2; }
    }
    return 1.45;
  }
  /** Drop test start: legs hanging, lowest tyre h metres above flat ground at y = 0. */
  function initialDropState(V, h, ph, vx) {
    const st = new Float64Array(14), c = Math.cos(ph), s = Math.sin(ph); st[0] = 10; st[2] = ph; let lowest = 1e9;
    for (let i = 0; i < 2; i++) {
      const lg = V.legs[i], th = airborneLegAngle(V, lg, ph); st[3 + i] = th;
      const K = legKinematics(lg, th); lowest = Math.min(lowest, s * K.bx + c * K.by);
    }
    st[1] = h + V.r - lowest; st[5] = vx; st[10] = st[11] = vx / V.r;
    return st;
  }
  /** Why a terrain run should end ('' = keep going). ctl = {u, stall} carries the stall timer. */
  function runStopReason(sim, V, lenM, ctl, dt) {
    const st = sim.st, xs = wheelXs(V, st);
    if (sim.status === 'diverged') return 'diverged';
    if (Math.max(xs[0], xs[1]) > lenM - V.r - 0.02) return 'finished';
    if (Math.min(xs[0], xs[1]) < V.r + 0.02 && st[5] < 0) return 'reached the start of the map';
    if (Math.abs(st[2]) > 1.2) return 'tipped over';
    if (Math.abs(ctl.u) > 0.05 && Math.abs(st[5]) < 0.01) ctl.stall += dt; else ctl.stall = 0;
    if (ctl.stall > 3) return 'stalled';
    return '';
  }


  /* =========================================================================
   * 6. TUNING  (spring rate, preload and damping that minimise CG acceleration)
   * =========================================================================
   * Candidate p = { k [N/mm], dh [mm], cb [N·s/m], cr [N·s/m] }, the same at both ends.
   * Preload is expressed as the ride height it produces: dh = static wheel travel
   * from the design start angle (+ = the rover sits lower, more sag). For each leg
   *     θ' = asin(sinθ0 − dh/L)
   * and buildVehicle solves the spring force that holds the rover level at θ'
   * (F = IR·(FN − (m_w + m_l/2)·g)), which fixes the preload x0 = F/k. Each dh
   * therefore maps to exactly one preload per leg, and the body stays level at rest.
   * The collar preload a mechanic sets (spring compression with the shock fully
   * extended) is x0 + ℓ_ride − ℓ_max.
   */
  const now = typeof performance !== 'undefined' && performance.now ? () => performance.now() : () => Date.now();

  /** Design variant for candidate p. Travel limits shift with the ride height so
   *  the absolute travel window relative to the body is unchanged. */
  function tunedDesign(design, p) {
    const d = JSON.parse(JSON.stringify(design));
    for (const s of ['L', 'R']) {
      const leg = d.legs[s];
      leg.th = Math.asin(clamp(Math.sin(leg.th * D2R) - p.dh / leg.L, -0.999, 0.999)) * R2D;
      leg.hmax -= p.dh; leg.droop += p.dh;
    }
    d.mode = 'manual'; d.kman = p.k;
    d.dyn = Object.assign({}, d.dyn, { cb: p.cb, cr: p.cr });
    return d;
  }
  /** Vehicle for candidate p plus the preload numbers to report (mm, N). */
  function tunedVehicle(design, p) {
    const d = tunedDesign(design, p), M = staticModel(d, { L: 0, R: 0 });
    if (!M.ok) return { err: 'That ride height puts a leg outside its reachable travel.' };
    const V = buildVehicle(d, M);
    if (V.err) return { err: V.err };
    const legs = V.legs.map(lg => ({
      th: lg.th0 * R2D, rideLen: lg.l0 * 1000, x0: lg.x0 * 1000,
      collar: (lg.x0 + lg.l0 - lg.lmax) * 1000, collarF: lg.k * (lg.x0 + lg.l0 - lg.lmax), springAtRide: lg.F0,
      inside: lg.l0 >= lg.lmin && lg.l0 <= lg.lmax, MR0: lg.MR0
    }));
    return { V, design: d, info: { k: p.k, dh: p.dh, cb: p.cb, cr: p.cr, legs } };
  }

  /**
   * Step-wise scenario runner (so a page can run it in small time slices).
   * kind 'run':  cfg = { T: contact terrain, lenM, x0 [m], v0 [m/s], u, tmax [s] }
   * kind 'drop': cfg = { h [m], ph [rad], vx [m/s], t [s] }
   * advance(ms) integrates for up to ms of wall time; returns true when finished.
   */
  function scenarioRunner(kind, V, cfg) {
    const run = kind === 'run';
    const sim = run ? new Simulation(V, cfg.T, initialRunState(V, cfg.T, cfg.x0, cfg.v0))
                    : new Simulation(V, flatTerrain(), initialDropState(V, cfg.h, cfg.ph, cfg.vx), { recDt: 0.001 });
    const tEnd = run ? cfg.tmax : cfg.t, u = run ? cfg.u : 0, ctl = { u, stall: 0 }, chunk = run ? 0.02 : 0.01;
    return {
      sim, kind, done: false, reason: '',
      progress() { return Math.min(1, sim.t / tEnd); },
      advance(ms) {
        const t0 = now();
        while (!this.done && now() - t0 < ms) {
          const tn = Math.min(tEnd, sim.t + chunk);
          if (!sim.advance(tn, u)) { this.reason = 'diverged'; this.done = true; break; }
          if (run) { const r = runStopReason(sim, V, cfg.lenM, ctl, chunk); if (r) { this.reason = r; this.done = true; break; } }
          if (sim.t >= tEnd - 1e-9) { this.reason = run ? 'reached the time limit' : 'done'; this.done = true; }
        }
        if (this.done && sim.rec.t[sim.rec.t.length - 1] < sim.t - 1e-9) sim.record(u);
        return this.done;
      }
    };
  }

  /**
   * CG acceleration metrics from a recording (g units, interval-averaged channels).
   *   a(t) = nz − 1   (vertical; 0 at rest), or |(nz − 1, nx)| with foreAft
   *   rms  = √(∫a² dt / T)
   *   peak = max |a| after a 10 ms moving average
   * fromContact: start at the first tyre contact (drop test, so free fall is not counted).
   * Also returns the closest approach to the bump stops (minMargin, metres; < 0 = bottomed).
   */
  function accelMetrics(R, V, opt) {
    opt = opt || {};
    const n = R.t.length;
    let i0 = 1;
    if (opt.fromContact) { i0 = R.N0.findIndex((v, i) => v > 0 || R.N1[i] > 0); if (i0 < 1) i0 = i0 < 0 ? n : 1; }
    let sum2 = 0, T = 0, peak = 0, peakT = NaN;
    const win = [], K = Math.max(1, Math.round(0.01 / Math.max(1e-6, (R.t[n - 1] - R.t[0]) / Math.max(1, n - 1))));
    let ws = 0;
    for (let i = i0; i < n; i++) {
      const dt = R.t[i] - R.t[i - 1], az = R.nza[i] - 1;
      const a = opt.foreAft ? Math.hypot(az, R.nxa[i]) : az;
      sum2 += a * a * dt; T += dt;
      win.push(a); ws += a; if (win.length > K) ws -= win.shift();
      if (win.length === K || i === n - 1) { const m = Math.abs(ws / win.length); if (m > peak) { peak = m; peakT = R.t[i]; } }
    }
    let minMargin = Infinity;
    for (let i = 0; i < n; i++) minMargin = Math.min(minMargin, R.l0[i] - V.legs[0].lmin, R.l1[i] - V.legs[1].lmin);
    return { rms: T > 0 ? Math.sqrt(sum2 / T) : 0, peak, peakT, minMargin, duration: T };
  }

  /**
   * Tuner objective (dimensionless, current setup = 1 per scenario).
   *   per scenario:  metric 'peak' → peak/peak₀, 'rms' → rms/rms₀, 'both' → mean of the two
   *   + 1 + shortfall/10 mm   if the coilovers come closer than minClear to the bump stops
   *   + 10                if the run fails where the current setup did not (stall, tip-over, unstable)
   *   total = w·J_run + (1 − w)·J_drop   (or the single scenario)
   */
  function tunerScore(res, base, cfg) {
    let J = 0;
    for (const s of ['run', 'drop']) {
      const m = res[s]; if (!m) continue;
      const b = base[s], w = res.run && res.drop ? (s === 'run' ? cfg.w : 1 - cfg.w) : 1;
      const nr = m.rms / Math.max(b.rms, 0.005), np = m.peak / Math.max(b.peak, 0.01);
      let js = cfg.metric === 'peak' ? np : cfg.metric === 'rms' ? nr : (np + nr) / 2;
      const mc = cfg.minClear != null ? cfg.minClear : 0;
      if (m.minMargin < mc) js += 1 + (mc - m.minMargin) * 1000 / 10;
      if (m.reason === 'diverged' || (s === 'run' && m.reason !== b.reason && m.reason !== 'finished')) js += 10;
      J += w * js;
    }
    return J;
  }

  function mulberry32(a) { return function () { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }

  /**
   * Bounded derivative-free search in normalized coordinates x ∈ [0,1]^dim, as a
   * generator: it yields points to evaluate and receives their objective values.
   *   1. the start point x0 (current setup)
   *   2. a Latin-hypercube scan of max(6, 2·dim + 2) points
   *   3. Nelder–Mead from the best point (reflect 1, expand 2, contract 0.5,
   *      shrink 0.5, initial simplex 0.15), candidates clamped to the box
   * Stops at the evaluation budget or when the simplex is smaller than 0.004.
   * Returns { x, f } of the best point.
   */
  function* tunerSearch(dim, x0, budget, seed) {
    const rng = mulberry32(seed || 12345), box = x => x.map(v => clamp(v, 0, 1));
    let n = 0, best = null;
    function* ev(x) { const f = yield x.slice(); n++; if (!best || f < best.f) best = { x: x.slice(), f }; return f; }
    yield* ev(box(x0));
    if (n >= budget) return best;
    const nInit = Math.min(budget - n, Math.max(6, 2 * dim + 2));
    const perms = [];
    for (let d = 0; d < dim; d++) { const p = Array.from({ length: nInit }, (_, i) => i); for (let i = nInit - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [p[i], p[j]] = [p[j], p[i]]; } perms.push(p); }
    for (let j = 0; j < nInit; j++) { yield* ev(perms.map(p => (p[j] + rng()) / nInit)); if (n >= budget) return best; }
    let S = [{ x: best.x.slice(), f: best.f }];
    for (let d = 0; d < dim; d++) {
      const v = best.x.slice(); v[d] = clamp(v[d] + (v[d] > 0.5 ? -0.15 : 0.15), 0, 1);
      S.push({ x: v, f: yield* ev(v) }); if (n >= budget) return best;
    }
    const comb = (a, b, t) => box(a.map((v, i) => v + t * (b[i] - v)));
    while (n < budget) {
      S.sort((a, b) => a.f - b.f);
      let diam = 0; for (const p of S) diam = Math.max(diam, Math.max(...p.x.map((v, i) => Math.abs(v - S[0].x[i]))));
      if (diam < 0.004) break;
      const w = S[S.length - 1], c = new Array(dim).fill(0);
      for (let i = 0; i < S.length - 1; i++) for (let d = 0; d < dim; d++) c[d] += S[i].x[d] / (S.length - 1);
      const xr = comb(c, w.x, -1), fr = yield* ev(xr); if (n >= budget) break;
      if (fr < S[0].f) {
        const xe = comb(c, w.x, -2), fe = yield* ev(xe);
        S[S.length - 1] = fe < fr ? { x: xe, f: fe } : { x: xr, f: fr };
      } else if (fr < S[S.length - 2].f) {
        S[S.length - 1] = { x: xr, f: fr };
      } else {
        const outside = fr < w.f, xc = outside ? comb(c, xr, 0.5) : comb(c, w.x, 0.5), fc = yield* ev(xc);
        if (fc < Math.min(fr, w.f)) S[S.length - 1] = { x: xc, f: fc };
        else {
          for (let i = 1; i < S.length && n < budget; i++) { const xs = comb(S[0].x, S[i].x, 0.5); S[i] = { x: xs, f: yield* ev(xs) }; }
        }
      }
    }
    return best;
  }

  return {
    G, D2R, R2D,
    // statics
    legPose, levelWheelLoads, springForceStatic, travelForLength, staticModel, travelSweep,
    // terrain
    smoothstep, OBSTACLES, presetMap, monotoneSpline, buildTerrain, contactTerrain, flatTerrain,
    // vehicle + dynamics
    RIGID_KT, buildVehicle, legKinematics, coiloverForce, solveLinear, prepareVehicle, assemble, derivatives,
    massMatrixAndForces, totalEnergy, newInfo, RECORD_KEYS, Simulation,
    // scenarios
    wheelXs, initialRunState, airborneLegAngle, initialDropState, runStopReason,
    // tuning
    tunedDesign, tunedVehicle, scenarioRunner, accelMetrics, tunerScore, tunerSearch
  };
});

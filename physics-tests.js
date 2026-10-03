/*
 * Rover Suspension Sizer — physics checks
 * ========================================
 * Independent checks of the math in physics.js. Each check compares the model
 * against something it was not built from: a force/moment balance, a finite
 * difference, an energy budget or a closed-form result.
 *
 *   Node:     node physics-tests.js
 *   Browser:  open verify.html (same folder as physics.js)
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory;
  else root.runPhysicsTests = factory;
})(typeof self !== 'undefined' ? self : this, function (RP) {
  'use strict';
  const D2R = Math.PI / 180;
  const clone = o => JSON.parse(JSON.stringify(o));
  const cross = (a, b) => a[0] * b[1] - a[1] * b[0];
  const results = [];
  const check = (group, name, value, limit, detail, cmp) => {
    const pass = cmp ? cmp(value, limit) : Math.abs(value) <= limit;
    results.push({ group, name, value, limit, pass, detail: detail || '' });
  };

  // Default design (same as the page's defaults)
  const LEG = { L: 600, th: 20, x1: 260, y1: 190, a: 305, e: 85, hmax: 100, droop: 50, ext: 255, stroke: 90 };
  const DESIGN = {
    m: 50, share: 50, cgx: 0, cgy: 0, rw: 120, mode: 'droop', G: 2.5, f: 2.0, kman: 8,
    legs: { L: { ...LEG }, R: { ...LEG } },
    dyn: { mw: 2.5, mleg: 1.0, ky: 250, J: 0.03, tyre: 'rigid', kt: 25, ct: 150, mu: 0.8, crr: 0.06, cb: 400, cr: 800, kbs: 300, Ts: 12, rpm: 150, dFront: true, dRear: true, coast: false }
  };
  // A deliberately lopsided design for the statics checks
  const ASYM = clone(DESIGN);
  Object.assign(ASYM, { m: 80, cgx: 40 });
  Object.assign(ASYM.legs.L, { L: 650, th: 25, x1: 180, y1: 170, a: 280, e: 40 });
  Object.assign(ASYM.legs.R, { x1: 200, y1: 150, a: 250, e: 0 });

  /* ---------------- 1. Statics ---------------- */
  for (const [label, des, trav] of [['symmetric, ride', DESIGN, { L: 0, R: 0 }], ['lopsided, front +60 mm', ASYM, { L: 10, R: 60 }]]) {
    const M = RP.staticModel(des, trav);
    if (!M.ok) { check('Statics', 'model solves (' + label + ')', 1, 0, M.msgs.join(' ')); continue; }
    let worstF = 0, worstM = 0;
    // body free-body: weight at CG, spring reaction at mount 1, pin reaction at O (global frame, mm & N)
    let BF = [0, -M.Wt], BM = cross([des.cgx, 0], [0, -M.Wt]);
    for (const s of ['L', 'R']) {
      const r = M.r[s], q = r.gn, sx = s === 'L' ? -1 : 1, Fs = r.Fh;
      // leg free-body in the outward frame
      const sumF = [Fs * q.u[0] + r.pin[0], r.FN + Fs * q.u[1] + r.pin[1]];
      const sumM = cross(q.W, [0, r.FN]) + cross(q.P2, [Fs * q.u[0], Fs * q.u[1]]);
      worstF = Math.max(worstF, Math.hypot(sumF[0], sumF[1])); worstM = Math.max(worstM, Math.abs(sumM));
      const P1 = [sx * des.legs[s].x1, des.legs[s].y1];
      const fSpringOnBody = [-Fs * q.u[0] * sx, -Fs * q.u[1]];
      const fPinOnBody = [-r.pin[0] * sx, -r.pin[1]];
      BF = [BF[0] + fSpringOnBody[0] + fPinOnBody[0], BF[1] + fSpringOnBody[1] + fPinOnBody[1]];
      BM += cross(P1, fSpringOnBody);
    }
    check('Statics', 'Leg ΣF = 0 (' + label + ')', worstF, 1e-9, 'N, worst leg');
    check('Statics', 'Leg ΣM about pivot = 0 (' + label + ')', worstM, 1e-6, 'N·mm, worst leg');
    check('Statics', 'Body ΣF = 0 (' + label + ')', Math.hypot(BF[0], BF[1]), 1e-9, 'N; weight + spring + pin reactions');
    check('Statics', 'Body ΣM about pivot = 0 (' + label + ')', BM, 1e-6, 'N·mm');
    // motion ratio: finite difference of spring length vs 1/leverage from the force balance
    let worst = 0;
    for (const s of ['L', 'R']) {
      const p = des.legs[s];
      for (const h of [-30, 0, 45, 90]) {
        const e = 1e-4, fd = -(RP.legPose(p, h + e).l - RP.legPose(p, h - e).l) / (2 * e);
        worst = Math.max(worst, Math.abs(fd - 1 / RP.legPose(p, h).IR));
      }
    }
    check('Statics', 'Motion ratio dΔ/dh = 1/leverage (' + label + ')', worst, 1e-6, 'virtual work, worst of 8 poses');
  }
  {
    const d1 = clone(DESIGN); d1.mode = 'droop';
    const M1 = RP.staticModel(d1, { L: 0, R: 0 });
    check('Sizing', 'Droop sizing: spring force at full droop = 0', M1.r.R.Fd, 1e-6, 'N');
    const d2 = clone(DESIGN); d2.mode = 'G'; d2.G = 2.7;
    const M2 = RP.staticModel(d2, { L: 0, R: 0 });
    check('Sizing', 'G sizing: wheel load at max bump = G × static', M2.r.R.Fb / M2.r.R.gb.IR / M2.r.R.FN0 - 2.7, 1e-9, '');
    const d3 = clone(DESIGN); d3.mode = 'freq'; d3.f = 1.8;
    const M3 = RP.staticModel(d3, { L: 0, R: 0 });
    check('Sizing', 'Frequency sizing: ride frequency = target', M3.r.R.f - 1.8, 1e-9, 'Hz');
  }

  /* ---------------- 2. Terrain contact ---------------- */
  {
    const flat = RP.contactTerrain(new Float64Array(401), 0.005), out = new Float64Array(9);
    const n1 = flat.contacts(1.0, 0.115, 0.12, out);
    check('Contact', 'Flat ground: one contact patch', n1, 1, 'patches', (v, l) => v === l);
    check('Contact', 'Flat ground: contact directly below the axle', Math.abs(out[0] - 1.0) + Math.abs(out[2] - 0.115), 1e-9, 'm');
    const step = RP.buildTerrain({ pts: [[0, 0], [2000, 0]], obs: [{ type: 'stepup', x: 1000, h: 80, w: 10, n: 1 }] }).ter;
    const n2 = step.contacts(0.9, 0.118, 0.12, out);
    check('Contact', 'Wheel against a step: floor + edge = two patches', n2, 2, 'patches', (v, l) => v === l);
  }

  /* ---------------- 3. Dynamics ---------------- */
  const M0 = RP.staticModel(DESIGN, { L: 0, R: 0 });
  const mkV = (patch, dynPatch) => { const d = clone(DESIGN); Object.assign(d, patch || {}); Object.assign(d.dyn, dynPatch || {}); return RP.buildVehicle(d, RP.staticModel(d, { L: 0, R: 0 })); };
  const air = RP.contactTerrain(new Float64Array(3001), 0.01);   // ground at y = 0; states below are 2 m up

  // Random states in the air
  let seed = 7; const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
  const randomState = () => { const st = new Float64Array(14); st[0] = 5; st[1] = 2 + rnd(); st[2] = (rnd() - 0.5) * 0.6; st[3] = 0.2 + rnd() * 0.4; st[4] = 0.2 + rnd() * 0.4; return st; };

  {
    const V = mkV({ cgx: 40, cgy: 90 }, { coast: true }); RP.prepareVehicle(V);
    let asym = 0, minPivot = Infinity, keErr = 0;
    for (let k = 0; k < 20; k++) {
      const st = randomState();
      for (let i = 5; i < 12; i++) st[i] = (rnd() - 0.5) * 4;
      const { M } = RP.massMatrixAndForces(V, air, st, 0);
      for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++) asym = Math.max(asym, Math.abs(M[i * 5 + j] - M[j * 5 + i]));
      // Cholesky to check positive definiteness
      const L = new Float64Array(25);
      for (let i = 0; i < 5; i++) for (let j = 0; j <= i; j++) {
        let s = M[i * 5 + j]; for (let p = 0; p < j; p++) s -= L[i * 5 + p] * L[j * 5 + p];
        if (i === j) { minPivot = Math.min(minPivot, s); L[i * 5 + i] = Math.sqrt(Math.max(s, 1e-300)); } else L[i * 5 + j] = s / L[j * 5 + j];
      }
      // ½ q̇ᵀ M q̇ + ½ J Σω²  vs  kinetic energy summed point by point
      let KEm = 0; for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++) KEm += 0.5 * st[5 + i] * M[i * 5 + j] * st[5 + j];
      KEm += 0.5 * V.J * (st[10] ** 2 + st[11] ** 2);
      keErr = Math.max(keErr, Math.abs(KEm - RP.totalEnergy(V, air, st).KE) / RP.totalEnergy(V, air, st).KE);
    }
    check('Equations of motion', 'Mass matrix symmetric', asym, 1e-12, 'kg·m², worst entry, 20 random states');
    check('Equations of motion', 'Mass matrix positive definite', minPivot, 0, 'smallest Cholesky pivot', (v) => v > 0);
    check('Equations of motion', '½q̇ᵀMq̇ equals the kinetic energy of every mass', keErr, 1e-12, 'relative error, 20 random states');

    // Conservative generalized forces = −∂PE/∂q (finite differences): gravity, springs, stops
    let worstQ = 0;
    for (let k = 0; k < 20; k++) {
      const st = randomState(); const { Q } = RP.massMatrixAndForces(V, air, st, 0);
      for (let i = 0; i < 5; i++) {
        const e = 1e-6, a = Float64Array.from(st), b = Float64Array.from(st); a[i] += e; b[i] -= e;
        const dPE = (RP.totalEnergy(V, air, a).PE - RP.totalEnergy(V, air, b).PE) / (2 * e);
        worstQ = Math.max(worstQ, Math.abs(Q[i] + dPE));
      }
    }
    check('Equations of motion', 'Generalized forces Q = −∂PE/∂q', worstQ, 1e-4, 'N or N·m, worst of 100 components (finite difference)');
  }

  // Free flight, no dissipation: energy must stay constant and converge with step size
  {
    const V = mkV({ cgx: 50, cgy: 100 }, { cb: 0, cr: 0, mu: 0, crr: 0, coast: true });
    V.cbs = 0;
    const st = randomState(); st[5] = 0.3; st[6] = 1; st[7] = 1.5; st[8] = 0.5; st[9] = -0.4; st[1] = 20;
    const spread = dtMax => {
      const S = new RP.Simulation(clone(V), air, st, { dtMax }); const E0 = RP.totalEnergy(S.V, air, S.st).E; let lo = E0, hi = E0;
      for (let k = 1; k <= 20; k++) { S.advance(k * 0.025, 0); const E = RP.totalEnergy(S.V, air, S.st).E; lo = Math.min(lo, E); hi = Math.max(hi, E); }
      return (hi - lo) / Math.abs(E0);
    };
    const s1 = spread(4e-4), s2 = spread(1e-4);
    check('Energy', 'Free flight, legs swinging: energy drift (dt 0.1 ms)', s2, 1e-4, 'relative, 0.5 s');
    check('Energy', 'Energy drift shrinks with smaller steps', s2 / s1, 1, 'ratio dt 0.1 ms / dt 0.4 ms', (v, l) => v < l);
  }
  // Drop onto the ground with every damper, friction and rolling loss switched off
  {
    const V = mkV({ cgx: 50, cgy: 100 }, { tyre: 'compliant', ct: 0, cb: 0, cr: 0, mu: 0, crr: 0, coast: true });
    V.cbs = 0;
    const S = new RP.Simulation(V, air, RP.initialDropState(V, 0.2, 3 * D2R, 0), { dtMax: 1e-4 });
    const E0 = RP.totalEnergy(V, air, S.st).E; let lo = E0, hi = E0, minl = Infinity;
    for (let k = 1; k <= 40; k++) { S.advance(k * 0.025, 0); const E = RP.totalEnergy(V, air, S.st).E; lo = Math.min(lo, E); hi = Math.max(hi, E); minl = Math.min(minl, ...S.rec.l0); }
    check('Energy', 'Lossless drop (tyre + springs + bump stops): energy drift', (hi - lo) / Math.abs(E0), 5e-4, 'relative, 1 s, min coilover ' + (minl * 1000).toFixed(1) + ' mm');
  }
  // Drop with the default (rigid) tyre: ground impulse = weight × time + momentum change
  {
    const V = mkV({}, {});
    const S = new RP.Simulation(V, air, RP.initialDropState(V, 0.3, 0, 0), { recDt: 0.001 });
    S.advance(2.0, 0); S.record(0);
    const R = S.rec, n = R.t.length; let J = 0;
    for (let i = 1; i < n; i++) J += (R.N0a[i] + R.N1a[i]) * (R.t[i] - R.t[i - 1]);
    // vertical momentum of every mass at the end
    const st = S.st, c = Math.cos(st[2]), sn = Math.sin(st[2]);
    const rcx = c * V.cg[0] - sn * V.cg[1];
    let py = V.mb * (st[6] + rcx * st[7]);
    for (let i = 0; i < 2; i++) { const K = RP.legKinematics(V.legs[i], st[3 + i]); const rbx = c * K.bx - sn * K.by, rdy = sn * K.dbx + c * K.dby;
      const vay = st[6] + rbx * st[7] + rdy * st[8 + i]; py += V.mw * vay + V.ml * (st[6] + (rbx * st[7] + rdy * st[8 + i]) / 2); }
    const m = V.mb + 2 * (V.mw + V.ml), expect = m * V.g * R.t[n - 1] + py;
    check('Momentum', 'Rigid-tyre drop: ∫N dt = m·g·t + final momentum', (J - expect) / expect, 1e-6, 'relative; impulse ' + J.toFixed(2) + ' N·s, expected ' + expect.toFixed(2));
  }
  // Rest on flat ground (rigid tyre): no drift, wheel loads = static
  {
    const V = mkV({}, { coast: true });
    const st = RP.initialRunState(V, air, 5, 0), y0 = st[1];
    const S = new RP.Simulation(V, air, st); S.advance(3, 0);
    check('Equilibrium', 'Rover at rest stays put (3 s)', Math.abs(S.st[1] - y0) * 1000 + Math.abs(S.st[0] - 5) * 1000, 0.01, 'mm drift');
    check('Equilibrium', 'Legs stay at the start angle', Math.abs(S.st[3] - V.legs[0].th0) + Math.abs(S.st[4] - V.legs[1].th0), 1e-5, 'rad');
    check('Equilibrium', 'Wheel load at rest = static load', Math.abs(S.rec.N1.at(-1) - V.FN[1]), 0.05, 'N');
  }
  // Steady speed on flat: motor torque = rolling-resistance moment
  {
    const V = mkV({}, {});
    const flat = RP.contactTerrain(new Float64Array(8001), 0.005);
    const S = new RP.Simulation(V, flat, RP.initialRunState(V, flat, 2, 0)); S.advance(12, 1);
    const N = V.FN[1], Mrr = V.crr * N * V.r, rl = V.r - N / V.kt;
    const vExp = V.w0 * rl * (1 - Mrr / V.Ts);
    check('Drive', 'Top speed on flat = motor curve vs rolling resistance', (S.st[5] - vExp) / vExp, 0.002, 'relative; expected ' + vExp.toFixed(4) + ' m/s, got ' + S.st[5].toFixed(4));
  }
  /* ---------------- 4. Tuner ---------------- */
  {
    const tv = RP.tunedVehicle(DESIGN, { k: 6, dh: 10, cb: 400, cr: 800 }), V = tv.V;
    const st = RP.initialRunState(V, air, 5, 0), S = new RP.Simulation(V, air, st); S.advance(2, 0);
    const hF = DESIGN.legs.R.L * (Math.sin(DESIGN.legs.R.th * D2R) - Math.sin(S.st[4]));
    check('Tuner', 'Candidate with +10 mm ride height rests there', Math.abs(hF - 10), 1e-6, 'mm of wheel travel after 2 s');
    check('Tuner', 'Candidate rests without drift', Math.abs(S.st[1] - st[1]) * 1000, 1e-6, 'mm');
    const lg = V.legs[1], F = RP.coiloverForce(V, lg, lg.lmax, 0);
    check('Tuner', 'Collar preload: spring force at full extension = k × preload', Math.abs(F - Math.max(0, tv.info.legs[1].collarF)), 1e-9, 'N');
    // search on a known bowl: minimum at (0.3, 0.7, 0.55, 0.2)
    const xs = [0.3, 0.7, 0.55, 0.2], bowl = x => x.reduce((a, v, i) => a + (i + 1) * (v - xs[i]) ** 2, 0);
    const it = RP.tunerSearch(4, [0.9, 0.1, 0.9, 0.9], 150, 3); let r = it.next();
    while (!r.done) r = it.next(bowl(r.value));
    check('Tuner', 'Search finds the minimum of a 4-D bowl (150 trials)', Math.max(...r.value.x.map((v, i) => Math.abs(v - xs[i]))), 0.01, 'worst coordinate error');
    // metrics on synthetic recordings: constant 0.5 g and a 1 g sine
    const mk = f => { const R = { t: [], nza: [], nxa: [], N0: [], N1: [], l0: [], l1: [] }; for (let i = 0; i <= 2000; i++) { const t = i * 0.001; R.t.push(t); R.nza.push(1 + f(t)); R.nxa.push(0); R.N0.push(1); R.N1.push(1); R.l0.push(0.2); R.l1.push(0.2); } return R; };
    const Vm = { legs: [{ lmin: 0.165 }, { lmin: 0.165 }] };
    const mc = RP.accelMetrics(mk(() => 0.5), Vm), ms = RP.accelMetrics(mk(t => Math.sin(2 * Math.PI * 2 * t)), Vm);
    check('Tuner', 'Metrics: constant 0.5 g gives rms = peak = 0.5', Math.abs(mc.rms - 0.5) + Math.abs(mc.peak - 0.5), 1e-9, 'g');
    check('Tuner', 'Metrics: 1 g sine gives rms = 1/√2', Math.abs(ms.rms - Math.SQRT1_2), 2e-3, 'g');
  }
  return results;
});

if (typeof require !== 'undefined' && typeof module !== 'undefined' && require.main === module) {
  const RP = require('./physics.js');
  const res = module.exports(RP);
  let fail = 0;
  for (const r of res) {
    if (!r.pass) fail++;
    console.log((r.pass ? 'PASS ' : 'FAIL ') + r.group.padEnd(20) + r.name.padEnd(62) + ' value ' + Number(r.value).toExponential(3) + '  limit ' + r.limit + '  ' + r.detail);
  }
  console.log('\n' + (res.length - fail) + ' of ' + res.length + ' checks passed');
  process.exit(fail ? 1 : 0);
}

# Rover Suspension Sizer

Design tool for the rover's rocker-leg coilover suspension: 2D statics and spring sizing, terrain runs, drop tests and a spring/damper tuner.

**Open it:** https://mikostrzewa.github.io/rover-suspension-sim/
**Documentation:** https://mikostrzewa.github.io/rover-suspension-sim/docs.html
**Physics checks:** https://mikostrzewa.github.io/rover-suspension-sim/verify.html

## Files

| File | What it is |
|---|---|
| `rover-suspension-sizer.html` | The tool (interface only) |
| `physics.js` | All the math; also runs in Node |
| `physics-tests.js` | 34 independent checks (`node physics-tests.js`) |
| `verify.html` | Runs the checks in the browser |
| `docs.html` | Documentation: every parameter and the models behind it |
| `index.html` | Redirects the site address to the tool |

## Notes for the team

- Your designs, presets and screen layout are saved in **your own browser**. Share designs with **Export presets** / **Import** in the sidebar, or **Download JSON** on the Statics tab.
- The site updates automatically about a minute after a change is pushed to `main`.
- The default numbers were eyeballed from early CAD; replace them with measured values.

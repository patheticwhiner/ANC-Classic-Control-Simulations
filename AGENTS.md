# Codex Project Instructions

## Claude Skill Bridge

Project-specific working procedures live under `.claude/skills/`. Codex must reuse them proactively; do not wait for the user to invoke a skill explicitly.

Route work as follows:

- Simulation report generation or revision: read `.claude/skills/anc-simulation-report/SKILL.md`.
- ANC figures, PSDs, convergence plots, or report graphics: also read `.claude/skills/anc-plot/SKILL.md`.
- Report, README, artifact, or reproduction consistency review: also read `.claude/skills/rep-audit/SKILL.md`.
- Publication-readiness or reviewer-style evaluation: read `.claude/skills/academic-reviewer-framework/SKILL.md` and only the relevant reference files it names.

When multiple routes apply, use them in this order:

1. `anc-simulation-report` defines the experiment and report contract.
2. `anc-plot` defines the visual evidence contract.
3. `rep-audit` checks that claims, files, figures, and numeric artifacts agree.

Read only relevant `.claude/memory/` files. Treat memory as historical context, not as authoritative evidence when code or frozen outputs disagree.

## Report Completion Gate

A simulation report is incomplete unless every reported scenario states its purpose, setup, comparison, quantitative result, interpretation, and conclusion. Reports must distinguish calibration from held-out evaluation, expose failed or saturated cases, and explain information-structure differences between controllers.

Generated claims must come from current result structures, MAT files, or CSV files. Avoid hard-coded numeric conclusions in report generators when the value can be derived from frozen outputs.

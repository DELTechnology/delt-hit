# Working State — Nature Protocols Revision

Last updated: 2026-06-01

This file is the central status document for the manuscript revision. It consolidates the current state across:

- `paper/main.tex`
- `paper/content/*.tex`
- `paper/tasks.txt`
- `paper/notes.md`
- `paper/rebuttle-letter.md`
- `paper/revision.md`
- `paper/reviewes.rtf`

## Purpose

Use this file as the primary coordination document for the revision. The other files remain useful, but they currently mix:

- historical planning
- partially stale status notes
- manuscript-independent brainstorms
- response-letter drafting

This document tracks what is actually present in the manuscript now, what still needs work, and where the source documents disagree.

## Source Documents and Their Current Role

`paper/main.tex` and `paper/content/*.tex`
- Current manuscript source of truth.
- Most reliable place to verify whether a claimed change has actually been incorporated.

`paper/tasks.txt`
- Short final-pass checklist.
- Good for submission-readiness tasks, but too narrow to capture overall revision state by itself.

`paper/revision.md`
- Broad revision checklist and Nature Protocols formatting plan.
- Useful as a backlog, but several items marked complete there still need verification against the latest manuscript and outputs.

`paper/notes.md`
- High-level synthesis of reviewer/editor asks and project history.
- Helpful context, but some status statements reflect earlier planning stages rather than the manuscript as it stands today.

`paper/rebuttle-letter.md`
- Draft point-by-point response letter.
- Contains useful reviewer phrasing, but several responses are incomplete, vague, or stronger than what is safely supported by the current manuscript.

`paper/reviewes.rtf`
- Original editorial letter and reviewer comments.
- Primary source for what Nature Protocols actually asked for.

## Current Manuscript Snapshot

## 1. Overall structure

The manuscript already follows a mostly Nature Protocols-like top-level structure:

- Abstract
- Key points
- Introduction
- Experimental Setup
- Materials
- Procedure
- Troubleshooting
- Timing
- Anticipated results
- Data and code availability
- Appendix / Supplementary Information

This means the remaining work is not mainly about global section order. The biggest open formatting work is inside the Procedure and around cross-referencing, consistency, and submission assets.

## 2. Items clearly already incorporated in the manuscript

- Title punctuation has been removed.
- The abstract includes:
  - required skill level
  - protocol duration
  - GitHub link
- The standalone Technical Overview section appears to have been removed.
- The Procedure is already converted to a continuous step sequence (`1` to `21`).
- Timing and Troubleshooting sections exist as standalone sections.
- The manuscript now explicitly mentions dual-display support.
- The DELi comparison section is substantially expanded.
- The backmatter already cites a Zenodo DOI: `10.5281/zenodo.20447074`.
- The anticipated-results section already includes a `zscore/` output layout.

## 3. Items that still look incomplete or risky

- Procedure formatting is improved but likely still not fully compliant with the editor’s requested style:
  - phase headings are still present as Procedure subheadings
  - some optional flows are not yet clearly expressed as `A/B/C` options or boxed optional procedures
  - some explanatory prose still sits outside the tightest possible step/callout structure
- Materials still describe a Conda/pip installation workflow, which may need alignment with the actual supported install path and the latest release/distribution method.
- The response letter is not yet safely aligned with the manuscript and evidence.
- Figure, table, and supplementary-note cross-references still need a systematic verification pass.
- Output filenames and command examples need verification against the current code and supporting material.
- The z-score method is mentioned in the manuscript and outputs, but it still needs an explicit consistency check against the latest implementation and CLI behavior.

## Document Disagreements to Resolve

## 1. Zenodo / code citation status

There is a mismatch across files:

- `paper/content/08-backmatter.tex` already states that code and source data are archived on Zenodo with DOI `10.5281/zenodo.20447074`.
- `paper/rebuttle-letter.md` says the release will be created shortly.
- `paper/revision.md` and `paper/tasks.txt` still treat release / DOI update as open.

Working assumption:
- The DOI exists, but the manuscript, release, and wording should be verified together before submission.

## 2. Abstract timing status

There is a mismatch across files:

- `paper/revision.md` says the abstract timing sentence still contains a placeholder.
- `paper/content/00-frontmatter.tex` currently says: "the workflow typically requires 1 hour to several hours."

Working assumption:
- The placeholder has been replaced in the manuscript, but it may still need tightening if benchmark-derived numbers are available.

## 3. Re-analysis / benchmarking status

There is tension between:

- `paper/rebuttle-letter.md`, which claims that full analyses and supplementary tables are already done
- `paper/notes.md`, which presents some of these items as still pending or blocked
- `paper/revision.md`, which marks benchmark generation as complete but still leaves manuscript integration open

Working assumption:
- Much of the computational work has likely been run already, but the manuscript and response letter still need verification against the actual generated assets before those claims should be treated as final.

## 4. Procedure compliance status

There is also a mismatch between:

- `paper/revision.md`, which marks several formatting tasks complete
- the current `paper/content/04-procedure.tex`, which still looks only partially adapted to the editor’s exact style guidance

Working assumption:
- The Procedure has improved substantially, but should not yet be considered fully finalized.

## Open Work by Priority

## 1. Submission-critical verification

- Verify every supplementary note is cited correctly.
- Verify every table is cited in the text.
- Verify every figure is cited in the text.
- Verify that all manuscript output names match current code outputs.
- Verify that all mentions of z-score match the current implementation and naming.

These are the tasks already listed in `paper/tasks.txt` and should remain the immediate pre-submission checklist.

## 2. Response-letter alignment

- Reconcile every claim in `paper/rebuttle-letter.md` with:
  - the current manuscript text
  - actual supplementary assets
  - actual code/release status
- Avoid claiming items are "done" unless they are visible in the manuscript or supporting files.
- Replace vague responses such as "Done" or "Partially addressed" with specific manuscript changes.

## 3. Procedure and formatting cleanup

- Recheck whether optional branches should be converted into clearer boxed or lettered alternatives.
- Recheck whether phase headings in the Procedure are still acceptable.
- Ensure all troubleshooting-linked steps are tagged consistently.
- Ensure timing references map cleanly to the final step numbering.

## 4. Release / repository consistency

- Confirm whether the cited Zenodo DOI corresponds to the intended public release state.
- Confirm whether the GitHub installation instructions match the release strategy.
- Update manuscript wording if the project should now cite a tagged release rather than a moving GitHub target.

## Recommended Working Rules

Until submission, use the documents like this:

- `paper/working-state.md`: central status and decision log
- `paper/tasks.txt`: short execution checklist for final verification
- `paper/rebuttle-letter.md`: reviewer response draft only
- `paper/revision.md`: backlog of broader revision requirements
- `paper/main.tex` + `paper/content/*.tex`: canonical manuscript text

If a status claim in another file conflicts with the manuscript, trust the manuscript first and then update the coordination docs.

## Immediate Next Pass

The most efficient next pass is:

1. verify all cross-references in manuscript and supplementary notes
2. verify output filenames and z-score terminology against the current code/output structure
3. update `paper/rebuttle-letter.md` so every response matches the current manuscript exactly
4. prune or update stale items in `paper/revision.md` and `paper/tasks.txt`

That sequence should leave one clean manuscript state and one clean response-letter state instead of several competing partial trackers.

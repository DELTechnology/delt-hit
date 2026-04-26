# Chemical Properties API Should Match Enumeration API

## Problem

The library properties API should follow the same shape as the numeration API so users do not have to learn two different invocation models for closely related workflows.

Right now, chemical properties feels like a separate interface. That makes it harder to compose with the rest of the pipeline and less obvious how to run it from existing outputs.

## Proposed Change

Update the library properties API so it accepts:

1. A `counts.txt` input
2. A library name

This should mirror the numeration API as closely as possible, including argument naming and the general command structure.

## Expected Output

In addition to API alignment, the output layout should be standardized so the generated properties artifact is written to:

```text
properties/<name>.parquet
```

Where `<name>` is the provided library name.

## Why This Helps

- keeps the user-facing workflow consistent across numeration and chemical properties
- makes it easier to feed the result of counting directly into downstream property generation
- gives outputs a predictable location and naming scheme
- makes batch runs and automation simpler

## Recommended Direction

Implement chemical properties as a sibling interface to numeration:

1. Accept `counts.txt`
2. Accept a library name
3. Write the result to `properties/<name>.parquet`

If the current implementation already computes the right data, the main change may be API and output-path normalization rather than chemistry logic itself.

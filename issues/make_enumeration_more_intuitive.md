# Make Enumeration More Intuitive

## Problem

Library enumeration is currently flexible, but the authoring model is harder to understand than it needs to be.

Today, the synthesis logic is effectively split across two places:

1. The `B0` / `B1` / ... sheets, where each row can define chemistry through columns such as `educt`, `reaction`, and `product`
2. The `reaction_graph` sheet, which defines additional arbitrary reaction steps through `educt_1`, `educt_2`, `reaction`, and `product`

This is powerful, but from a user perspective it has a few drawbacks:

- Users must understand an internal graph model rather than only describing synthesis routes
- Additional steps are expressed at a very low level as graph edges
- Intermediate names such as `product_0`, `product_1B`, etc. become part of the authoring interface
- It is easy for the workbook to feel like it has multiple overlapping sources of truth
- The user has to mentally reconstruct which subsets of `B0`, `B1`, ... are compatible

The main issue is not enumeration itself. The main issue is that the current representation exposes the implementation model too directly.

## Important Constraint

Any simplification must preserve an important aspect of real libraries:

- different mutually exclusive subsets of `B0`, `B1`, ... often participate in different reactions
- we cannot assume that all members of `B0` use the same reaction step
- therefore a single uniform per-cycle topology is usually not expressive enough

This means a naive simplification like "all `B0` entries use the same step, all `B1` entries use the same step" will not work for many real libraries.

## Why Not Use Free-Form Route Strings

One possible alternative would be to define synthesis routes as strings, for example:

```text
synthesis_0: B0 + scaffold_0 -R1-> product_0 + B1 -R2-> final_0
synthesis_1: scaffold_pre -R9-> scaffold_0 + B0 -R2-> product_0 + B1 -R2-> final_1
```

This looks compact, but it likely creates more problems than it solves:

- users would have to write a mini language correctly
- validation would become harder
- explicit chemistry would be hidden inside strings
- the system would need reserved names like `B0`, `B1`, etc.
- parse errors and authoring mistakes would be harder to diagnose than in a table

For these reasons, free-form route strings are probably not the right direction.

## Potential Solution: Explicit Routes Plus Structured Steps

A more promising design is to introduce explicit route membership, while keeping route definitions structured in a table rather than in strings.

### High-level idea

- `B0`, `B1`, ... sheets define building blocks and the route(s) they belong to
- a new `reaction_steps` sheet defines the non-building-block chemistry and overall route topology
- the internal `reaction_graph` is no longer the primary authoring interface
- instead, the graph is compiled from the route definitions

### Example

`B0`

| codon | smiles | route |
|---|---|---|
| bb_a | ... | route_0 |
| bb_b | ... | route_0 |
| bb_c | ... | route_1 |

`B1`

| codon | smiles | route |
|---|---|---|
| bb_d | ... | route_0 |
| bb_e | ... | route_1 |

`reaction_steps`

| route | step | input_1 | input_2 | reaction | output |
|---|---:|---|---|---|---|
| route_0 | 1 | scaffold_0 | B0 | R1 | product_0 |
| route_0 | 2 | product_0 | B1 | R2 | final_0 |
| route_1 | 1 | scaffold_pre |  | R9 | scaffold_0 |
| route_1 | 2 | scaffold_0 | B0 | R7 | product_x |
| route_1 | 3 | product_x | B1 | R8 | final_1 |

In this model:

- `bb_a` and `bb_b` are enumerated with `route_0`
- `bb_c` is enumerated with `route_1`
- different routes can use different reactions
- route compatibility is explicit rather than inferred indirectly from intermediate names

## Why This Helps

This design keeps the full expressive power needed for real libraries, but it changes the user-facing abstraction:

- users think in terms of discrete synthesis routes
- route-specific additional steps are expressed as ordered reactions, not raw graph edges
- building-block sheets focus on building-block identity plus route membership
- the graph remains an internal execution model, not the main authoring surface

In other words:

- current model: "manually define a graph"
- proposed model: "define route members and route steps; derive the graph"

## Tradeoffs and Open Questions

This is not a free simplification. It introduces its own decisions:

- Should `route` allow one value or multiple values per building block row?
- How should shared steps across routes be represented?
- Are routes truly discrete for all relevant libraries, or are there cases where compatibility is more fluid?
- Do we still need some escape hatch equivalent to the current `reaction_graph` for unusual edge cases?

The last question is likely the most important. The current `reaction_graph` exists because arbitrary additional steps are genuinely hard to express otherwise. A route-based design may cover the common case much better, but there may still be rare cases where a lower-level mechanism is needed.

## Recommended Direction

Explore replacing `reaction_graph` as the primary user-facing model with:

1. A `route` column on `B0`, `B1`, ... sheets
2. A structured `reaction_steps` sheet containing route-specific ordered steps
3. Compilation of those steps into the current internal graph representation

This would likely make enumeration authoring more intuitive without giving up support for route-specific chemistry.

## Non-Goal

This proposal does **not** argue that the internal graph representation should be removed.

The goal is only to make the authoring model more intuitive:

- keep the graph internally if it remains the best execution model
- reduce the amount of graph bookkeeping users must do by hand

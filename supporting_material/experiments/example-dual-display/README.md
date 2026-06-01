# Example Dual-Display Configuration

This folder contains a minimal dual-display Excel configuration example:

- [example-dual-display.xlsx](/Users/adrianomartinelli/projects/delt-hit/supporting_material/experiments/example-dual-display/example-dual-display.xlsx)

Dual-display-specific columns:

- `structure` sheet: use the `strand` column on every `building_block` row and assign each row to either `A` or `B`.
- `structure` sheet: `reverse` and `complement` define how each DNA region is oriented in the sequencing construct.
- Building-block sheets (`B0`, `B1`, ...): keep `educt`, `reaction`, and `product` definitions within the same strand branch. DELT-Hit validates that an `A`-strand building block only connects to `A`-strand intermediates, and likewise for `B`.

Enumeration behavior:

- Single-display libraries produce one `smiles` column.
- Dual-display libraries are enumerated strand-wise and produce `smiles_a` and `smiles_b`.

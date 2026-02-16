# Phylogenetic Tree Data Provenance

## Source

All trees were downloaded from the **Open Tree of Life** synthetic tree via the v3 API on 2026-02-16.

- API endpoint: `https://api.opentreeoflife.org/v3/tree_of_life/subtree`
- Method: POST with `{"ott_id": <ID>}` for each clade
- Format: Newick (with OTT identifiers as internal node labels)

## Citation

Hinchliff, C. E., et al. (2015). Synthesis of phylogeny and taxonomy into a
comprehensive tree of life. *Proceedings of the National Academy of Sciences*,
112(41), 12764--12769. https://doi.org/10.1073/pnas.1423041112

Open Tree of Life: https://opentreeoflife.github.io/

## Inventory

| File                | Clade           | Common name                    | Rank       | Tips   | OTT ID  |
|---------------------|-----------------|--------------------------------|------------|--------|---------|
| homininae.nwk       | Homininae       | great apes and human ancestors | subfamily  |     24 |  312031 |
| delphinidae.nwk     | Delphinidae     | oceanic dolphins               | family     |     60 |  698406 |
| felidae.nwk         | Felidae         | cats                           | family     |     91 |  563159 |
| mustelidae.nwk      | Mustelidae      | weasels, otters, badgers       | family     |    106 |  348043 |
| equidae.nwk         | Equidae         | horses, zebras, asses          | family     |    171 |  698424 |
| salamandridae.nwk   | Salamandridae   | newts and salamanders          | family     |    202 |  566011 |
| canidae.nwk         | Canidae         | dogs, wolves, foxes            | family     |    211 |  770319 |
| accipitridae.nwk    | Accipitridae    | hawks, eagles, kites           | family     |    529 | 1036185 |
| primates.nwk        | Primates        | primates                       | order      |    748 |  913935 |
| cetacea.nwk         | Cetacea         | whales and dolphins            | order      |    848 |  622916 |
| papilionidae.nwk    | Papilionidae    | swallowtail butterflies        | family     |   1082 |  661439 |
| aves.nwk            | Aves            | birds                          | class      | 18,990 |   81461 |

## Notes

- Tip counts reflect leaf nodes (terminal taxa) in the synthetic tree. These include
  subspecies and some fossil taxa, so counts may exceed the number of extant species.
- The Homininae tree includes extinct hominin species (H. erectus, H. habilis, etc.).
- The Equidae tree includes many fossil taxa.
- The Aves tree is very large (18,990 tips, 733 KB). For computational experiments,
  consider using the Accipitridae subtree or sampling.
- Internal node labels in the Newick files carry OTT identifiers (e.g., `Panthera_ott563154`).
- Branch lengths are not included (topology only).

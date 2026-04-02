# Nematode 0.2.2 (2026-03-18)

## DATA UPDATES

-   **Updated `nematode.info`:** Refreshed taxonomic and functional traits dataset (2,484 --> 2524 genera) from Nemaplex.UCDavis.edu (Revision Date: 03/18/2026).

-   **Updated `nematode.ave.mass`:** Refreshed genus and family average body mass dataset (987 --> 1094 entries) from Nemaplex.UCDavis.edu (Revision Date: 03/18/2026).

## BUG FIXES

- Fixed various minor issues and improved overall stability.


# Nematode 0.2.1 (2025-09-07)

## BUG FIXES

-   **Formula correction in SRI calculation:**
    -   Fixed formula error in `cal.SRI()` affecting Species Richness Index (SRI) calculation
    -   Fixed related error in `Ecological.Indices()` function
    -   **Impact:** Incorrect formula produced erroneous SRI values
    -   Upgrade to latest version required for correct calculations

# Nematode 0.2.0 (2025-08-20)

## NEW FEATURES

-   Added function `cal.TD`: Trophic Diversity
-   Added function `cal.H`: Shannon-Wiener Index
-   Added function `cal.J`: Pielou's Evenness Index
-   Added function `cal.Simpson`: Simpson Index
-   Added function `cal.WI`: Wasilewska Index
-   Added function `cal.MI`: Maturity Index
-   Added function `cal.PPI`: Plant Parasite Index
-   Added function `cal.SRI`: Species Richness Index
-   Added function `cal.NCR`: Nematode Channel Ratio
-   Added function `cal.CI`: Channel Index
-   Added function `cal.BI`: Basic Index
-   Added function `cal.EI`: Enrichment Index
-   Added function `cal.SI`: Structure Index
-   Added function `Ecological.Indices`
-   Added function `NMF`
-   Added function `NEF`
-   Added function `check_nematode_genus`
-   Added function `diet_rel_abundance`
-   Added function `fuzzy_genus_match`
-   Added function `num_species`
-   Added function `rel_abundance`
-   Added function `runNMDS`
-   Added function `runPCoA`
-   Added function `runSimper`

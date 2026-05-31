# Plans

## Done
- **Fix `internal_node` clone assignment (tree-aware).** Root cause of the
  single-clone bug: `select_internal_nodes_by_branch_length` selected every
  internal node when `min_branch_length=0.0` (default), and `cut_tree_at_nodes`
  assigned in preorder so the root claimed all cells. Fixed by (1) excluding the
  root, (2) using a strict `>` threshold, (3) deriving an automatic threshold
  (mean + 1 SD of internal edge lengths) when `min_branch_length<=0`, and
  (4) assigning deepest clades first.

## Next
- **`distance` strategy with automatic threshold.** Implement auto-selection of
  `--distance_threshold` so users don't have to guess it. Approach: compute the
  pairwise tree distance matrix, build the Ward linkage, and pick the cut height
  automatically (e.g. largest linkage gap / elbow, or a silhouette/gap-statistic
  sweep over candidate cut heights). Fall back to the user-supplied
  `--distance_threshold` when provided. Goal: `--clone_strategy distance` works
  with no manual threshold and yields several biologically sensible clones.
- Consider an optional `--n_clones` target that both strategies can honor.

library(DiagrammeR)

grViz("
digraph brms_path {

  graph [layout = dot, rankdir = LR]
  node  [shape = box, fontsize = 11]

  light [label = 'lighttreatment']
  gcl   [label = 'log(gcl)']
  fg    [label = 'log(fgmax)']
  lam   [label = 'log(lambda_mean)']
  tau   [label = 'log(tau_mean)']

  node [shape = ellipse]
  acc  [label = 'accession (RE)']
  phy  [label = 'phylogeny (RE; A)']

  # Directed (regression) paths
  light -> gcl
  light -> fg
  light -> lam
  light -> tau
  gcl   -> lam
  gcl   -> tau
  fg    -> lam
  fg    -> tau

  # Random effects influence all responses
  acc -> gcl
  acc -> fg
  acc -> lam
  acc -> tau

  phy -> gcl
  phy -> fg
  phy -> lam
  phy -> tau

  # Residual correlations (from set_rescor(TRUE))
  edge [dir=both, arrowhead=none, arrowtail=none, style=dashed]
  lam -> tau
  lam -> gcl
  lam -> fg
  tau -> gcl
  tau -> fg
  gcl -> fg

  # Annotation about interactions
  labelloc = 't'
  label = 'Interactions: lighttreatment × log(gcl) × log(fgmax) for lambda and tau'
}
")

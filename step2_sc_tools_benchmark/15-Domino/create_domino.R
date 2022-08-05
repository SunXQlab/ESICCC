domino <- setClass(
  Class = 'domino',
  slots = c(
    db_info = 'list',
    z_scores = 'matrix',
    counts = 'AnyMatrix',
    clusters = 'factor',
    features = 'matrix',
    cor = 'matrix',
    linkages = 'list',
    clust_de = 'matrix',
    misc = 'list',
    cl_signaling_matrices = 'list',
    signaling = 'matrix'
  )
)

add_rl_column = function(map, map_ref, conv, new_name){
  map_in_ref = match(map[[map_ref]], conv[,1])
  not_in_ref = which(is.na(map_in_ref))
  if(length(not_in_ref > 0)){
    not_in_ref_map = cbind.data.frame(map[not_in_ref,], as.character(NA), 
                                      stringsAsFactors = FALSE)
    colnames(not_in_ref_map)[ncol(not_in_ref_map)] = new_name
    rownames(not_in_ref_map) = c()
  } else {
    not_in_ref_map = c()
  }
  new_map = c()
  for(r_id in 1:nrow(map)){
    row = map[r_id,]
    conv_ids = which(conv[,1] == row[[map_ref]])
    for(id in conv_ids){
      new_row = c(as.matrix(row), conv[id, 2])
      new_map = rbind(new_map, new_row)
    }
  }
  rownames(new_map) = c()
  colnames(new_map) = c(colnames(map), new_name)
  new_map = rbind.data.frame(new_map, not_in_ref_map, stringsAsFactors = FALSE)
  new_map = data.frame(new_map, stringsAsFactors = FALSE)
}

create_domino <- function (signaling_db, features, ser = NULL, counts = NULL, 
          z_scores = NULL, clusters = NULL, use_clusters = TRUE, df = NULL, 
          gene_conv = NULL, verbose = TRUE, use_complexes = TRUE, rec_min_thresh = 0.025, 
          remove_rec_dropout = TRUE, tf_selection_method = "clusters", 
          tf_variance_quantile = 0.5) 
{
  dom = domino()
  dom@misc[["tar_lr_cols"]] = c("R.orig", "L.orig")
  dom@misc[["create"]] = TRUE
  dom@misc[["build"]] = FALSE
  dom@misc[["build_vars"]] = NULL
  if (!is.null(ser) & (!is.null(clusters) | !is.null(z_scores) | 
                       !is.null(counts))) {
    warning("Ser and z_score, clusters, or counts provided. Defaulting to ser.")
  }
  if (is.null(ser) & (is.null(clusters) | is.null(z_scores) | 
                      is.null(counts))) {
    stop("Either ser or clusters and z_scores must be provided")
  }
  if (!(tf_selection_method %in% c("all", "clusters", "variable"))) {
    stop("tf_selection_method must be one of all, clusters, or variable")
  }
  if (verbose) {
    print("Reading in and processing signaling database")
  }
  genes = read.csv(paste0(signaling_db, "/genes.csv"), stringsAsFactors = FALSE)
  proteins = read.csv(paste0(signaling_db, "/proteins.csv"), 
                      stringsAsFactors = FALSE)
  interactions = read.csv(paste0(signaling_db, "/interactions.csv"), 
                          stringsAsFactors = FALSE)
  complexes = read.csv(paste0(signaling_db, "/complexes.csv"), 
                       stringsAsFactors = FALSE)
  dom@db_info = list(genes = genes, proteins = proteins, interactions = interactions, 
                     complexes = complexes)
  rec_uniprot = proteins$uniprot[which(proteins$receptor == 
                                         "True")]
  lig_uniprot = proteins$uniprot[which(proteins$receptor != 
                                         "True")]
  if (use_complexes) {
    rec_complexes = complexes$complex_name[which(complexes$receptor == 
                                                   TRUE)]
    lig_complexes = complexes$complex_name[which(complexes$receptor != 
                                                   TRUE)]
  }
  rl_map = matrix(0, ncol = 2, nrow = 0)
  colnames(rl_map) = c("R.uniprot", "L.uniprot")
  for (i in 1:nrow(interactions)) {
    a_pname = interactions$protein_name_a[i]
    b_pname = interactions$protein_name_b[i]
    if (!use_complexes) {
      if (a_pname == "" | b_pname == "") {
        next
      }
    }
    aid = interactions$partner_a[i]
    bid = interactions$partner_b[i]
    if (length(which(rec_uniprot == aid)) != 0) {
      recs = rec_uniprot[which(rec_uniprot == aid)]
    }
    else if (length(which(lig_uniprot == aid)) != 0) {
      ligs = lig_uniprot[which(lig_uniprot == aid)]
    }
    else if (length(which(rec_complexes == aid)) != 0) {
      comp = rec_complexes[which(rec_complexes == aid)]
      recs = complexes[which(complexes[, 1] == comp), 2:5]
      recs = recs[-which(recs == "")]
    }
    else if (length(which(lig_complexes == aid)) != 0) {
      comp = lig_complexes[which(lig_complexes == aid)]
      ligs = complexes[which(complexes[, 1] == comp), 2:5]
      ligs = ligs[-which(ligs == "")]
    }
    else {
      stop(paste("Partner A has no comp or prot match in row", 
                 i))
    }
    if (length(which(rec_uniprot == bid)) != 0) {
      recs = rec_uniprot[which(rec_uniprot == bid)]
    }
    else if (length(which(lig_uniprot == bid)) != 0) {
      ligs = lig_uniprot[which(lig_uniprot == bid)]
    }
    else if (length(which(rec_complexes == bid)) != 0) {
      comp = rec_complexes[which(rec_complexes == bid)]
      recs = complexes[which(complexes[, 1] == comp), 2:5]
      recs = recs[-which(recs == "")]
    }
    else if (length(which(lig_complexes == bid)) != 0) {
      comp = lig_complexes[which(lig_complexes == bid)]
      ligs = complexes[which(complexes[, 1] == comp), 2:5]
      ligs = ligs[-which(ligs == "")]
    }
    else {
      stop(paste("Partner B has no comp or prot match in row", 
                 i))
    }
    for (l in ligs) {
      for (r in recs) {
        rl_map = rbind(rl_map, c(r, l))
      }
    }
  }
  rec_prots = as.character(unique(rl_map[, 1]))
  rec_orig = genes[match(rec_prots, genes[, 2]), 3]
  rl_map = data.frame(rl_map)
  rl_map = add_rl_column(rl_map, "R.uniprot", cbind(rec_prots, 
                                                    rec_orig), "R.orig")
  lig_prots = as.character(unique(rl_map[, 2]))
  lig_orig = genes[match(lig_prots, genes[, 2]), 3]
  rl_map = add_rl_column(rl_map, "L.uniprot", cbind(lig_prots, 
                                                    lig_orig), "L.orig")
  if (!is.null(gene_conv)) {
    conv = convert_genes(rl_map$R.orig, from = gene_conv[1], 
                         to = gene_conv[2])
    rl_map = add_rl_column(rl_map, "R.orig", conv, "R.conv")
    conv = convert_genes(rl_map$L.orig, from = gene_conv[1], 
                         to = gene_conv[2])
    rl_map = add_rl_column(rl_map, "L.orig", conv, "L.conv")
    dom@misc[["tar_lr_cols"]] = c("R.conv", "L.conv")
  }
  rl_map = unique.data.frame(rl_map)
  rl_list = list()
  tar_map = rl_map[, dom@misc$tar_lr_cols]
  for (rec in unique(tar_map[, 1])) {
    if (is.na(rec)) {
      next
    }
    ligs = unique(tar_map[which(tar_map[, 1] == rec), 2])
    for (lig in ligs) {
      if (!is.na(lig)) {
        rl_list[[rec]] = c(rl_list[[rec]], lig)
      }
    }
  }
  dom@linkages[["rec_lig"]] = rl_list
  dom@misc[["rl_map"]] = rl_map
  if (verbose) {
    print("Getting z_scores, clusters, and counts")
  }
  if (!is.null(ser)) {
    z_scores = ser@assays$RNA@scale.data
    if (use_clusters) {
      clusters = ser@active.ident
    }
    counts = ser@assays$RNA@counts
  }
  dom@z_scores = z_scores
  if (!is.null(clusters)) {
    dom@clusters = clusters
  }
  if (class(features) == "character") {
    features = read.csv(features, row.names = 1, check.names = FALSE)
  }
  features = features[, colnames(dom@z_scores)]
  dom@features = as.matrix(features)
  if (tf_selection_method == "clusters") {
    p_vals = matrix(1, nrow = nrow(features), ncol = length(levels(dom@clusters)))
    rownames(p_vals) = rownames(features)
    colnames(p_vals) = levels(dom@clusters)
    if (verbose) {
      print("Calculating feature enrichment by cluster")
      clust_n = length(levels(dom@clusters))
    }
    for (clust in levels(dom@clusters)) {
      if (verbose) {
        cur = which(levels(dom@clusters) == clust)
        print(paste0(cur, " of ", clust_n))
      }
      cells = which(dom@clusters == clust)
      for (feat in rownames(dom@features)) {
        p_vals[feat, clust] = wilcox.test(dom@features[feat, 
                                                       cells], dom@features[feat, -cells], alternative = "g")$p.value
      }
    }
    dom@clust_de = p_vals
  }
  if (tf_selection_method == "all") {
    dom@clusters = factor()
  }
  if (tf_selection_method == "variable") {
    dom@clusters = factor()
    variances = apply(dom@features, 1, function(x) {
      sd(x)/mean(x)
    })
    keep_n = length(variances) * tf_variance_quantile
    keep_id = which(rank(variances) > keep_n)
    dom@features = dom@features[names(keep_id), ]
  }
  if (class(df) == "character") {
    df = read.csv(df, skip = 3, header = FALSE, stringsAsFactors = FALSE)
  }
  if (!is.null(df)) {
    tf_targets = list()
    for (row in 1:nrow(df)) {
      tf = df[row, 1]
      tf_genes = df[row, 9]
      tf_genes = gsub("[\\(']", "", regmatches(tf_genes, 
                                               gregexpr("\\('.*?'", tf_genes))[[1]])
      tf_targets[[tf]] = unique(c(tf_targets[[tf]], tf_genes))
    }
    dom@linkages[["tf_targets"]] = tf_targets
  }
  else {
    dom@linkages[["tf_targets"]] = NULL
  }
  dom@counts = counts
  all_receptors = unique(names(dom@linkages$rec_lig))
  zero_sum = Matrix::rowSums(counts == 0)
  keeps = which(zero_sum < 0.975 * ncol(counts))
  ser_receptors = intersect(names(keeps), all_receptors)
  rho = matrix(0, nrow = length(ser_receptors), ncol = nrow(dom@features))
  rownames(rho) = ser_receptors
  colnames(rho) = rownames(dom@features)
  if (verbose) {
    print("Calculating correlations")
    n_tf = nrow(dom@features)
  }
  for (module in rownames(dom@features)) {
    if (verbose) {
      cur = which(rownames(dom@features) == module)
      print(paste0(cur, " of ", n_tf))
    }
    if (!is.null(dom@linkages$tf_targets)) {
      tf = substring(module, first = 1, last = nchar(module) - 
                       3)
      module_targets = tf_targets[[tf]]
      module_rec_targets = intersect(module_targets, ser_receptors)
    }
    scores = dom@features[module, ]
    rhorow = rep(0, length(ser_receptors))
    names(rhorow) = ser_receptors
    for (rec in ser_receptors) {
      if (remove_rec_dropout) {
        keep_id = which(dom@counts[rec, ] > 0)
        rec_z_scores = dom@z_scores[rec, keep_id]
        tar_tf_scores = scores[keep_id]
      }
      else {
        rec_z_scores = dom@z_scores[rec, ]
        tar_tf_scores = scores
      }
      if (sum(tar_tf_scores) == 0) {
        rhorow[rec] = 0
        next
      }
      cor = cor.test(rec_z_scores, tar_tf_scores, method = "spearman", 
                     alternative = "greater")
      rhorow[rec] = cor$estimate
    }
    if (length(module_rec_targets > 0)) {
      rhorow[module_rec_targets] = 0
    }
    rho[, module] = rhorow
  }
  colnames(rho) = rownames(dom@features)
  dom@cor = rho
  return(dom)
}

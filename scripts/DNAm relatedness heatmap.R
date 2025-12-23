suppressPackageStartupMessages({
  library(openxlsx)
  library(tidyverse)
  library(lubridate)
  library(glue)
  library(kinship2)
  library(ggplot2)
  library(viridis)
  library(patchwork) 
})

# ---------------------------- CONFIG -----------------------------------------
IND_PATH        <- "./data/Peromyscus.xlsx"
MATING_PATH     <- "./data/Mating Records.xlsx"
DNAM_PATH       <- "./data/DNAm mice with EAA.xlsx"
TARGET_SPECIES  <- c("BW")

GENERATION_DEPTH <- 6

OUT_DIR         <- "./relatedness"
OUT_FILE_TMPL   <- "{species} - DNAm relatedness (gen {gen_depth}).xlsx"

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# --------------------------- HELPERS -----------------------------------------

clean_column_names <- function(df) {
  names(df) <- gsub("[^[:alnum:]_]", "", names(df))
  df
}

normalize_id <- function(x) {
  y <- gsub("[^[:alnum:]]", "", as.character(x))
  y <- ifelse(grepl("^0", y), paste0(substr(y, 2, nchar(y)), "00000"), y)
  suppressWarnings(as.numeric(y))
}

clean_mating <- function(x) {
  y <- gsub("[^[:alnum:]]", "", as.character(x))
  y <- toupper(y)
  ifelse(y == "", NA_character_, y)
}

normalize_stock <- function(x) toupper(trimws(gsub("[^[:alnum:] ]", "", as.character(x))))

pick_col <- function(df, candidates, required = TRUE, label = "") {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) {
    if (required) stop(glue("Missing required column(s) for {label}: {paste(candidates, collapse=", ")}"))
    return(NA_character_)
  }
  hit[[1]]
}

get_ancestors <- function(ids, ped_data, max_generations) {
  work <- unique(stats::na.omit(as.numeric(ids)))
  anc <- work
  for (i in seq_len(max_generations)) {
    if (length(work) == 0) break
    parents <- ped_data %>%
      filter(ID %in% work) %>%
      dplyr::select(Dam, Sire) %>%
      unlist(use.names = FALSE) %>%
      stats::na.omit() %>%
      unique()
    new_ids <- setdiff(parents, anc)
    if (length(new_ids) == 0) break
    anc  <- unique(c(anc, new_ids))
    work <- new_ids
  }
  anc
}

# --- FULLY CORRECTED FUNCTION ---
compute_parent_relatedness <- function(dam_id, sire_id, IND, DAMSIRE, ped_data, depth) {
  if (any(is.na(c(dam_id, sire_id)))) return(NA_real_)
  
  # Scope by ancestors up to depth
  base_ids <- unique(c(dam_id, sire_id))
  rel_ids  <- unique(c(base_ids, get_ancestors(base_ids, ped_data, depth)))
  
  # Start from IND rows for those IDs
  # Renamed 'rel' to 'rel_df'
  rel_df <- IND %>%
    filter(ID %in% rel_ids) %>%
    distinct(ID, .keep_all = TRUE) %>%
    left_join(DAMSIRE, by = "MatingNumber") %>%
    distinct(ID, .keep_all = TRUE) %>%
    mutate(
      Dam  = na_if(Dam, 0),
      Sire = na_if(Sire, 0),
      Sex_k2 = case_when(
        toupper(Sex) %in% c("M","MALE") ~ 1L,
        toupper(Sex) %in% c("F","FEMALE") ~ 2L,
        TRUE ~ 0L
      )
    )
  
  # Ensure the focal parents exist
  if (!(dam_id %in% rel_df$ID)) rel_df <- bind_rows(rel_df, tibble(ID = dam_id, MatingNumber = NA_character_, Sex = "F", Dam = NA_real_, Sire = NA_real_, Sex_k2 = 2L))
  if (!(sire_id %in% rel_df$ID)) rel_df <- bind_rows(rel_df, tibble(ID = sire_id, MatingNumber = NA_character_, Sex = "M", Dam = NA_real_, Sire = NA_real_, Sex_k2 = 1L))
  
  # Ensure *all* referenced parents exist as rows (founders if absent)
  all_parent_ids <- unique(stats::na.omit(c(rel_df$Dam, rel_df$Sire)))
  missing_parents <- setdiff(all_parent_ids, rel_df$ID)
  if (length(missing_parents) > 0) {
    # Infer sex where possible from being listed as Sire/Dam
    is_sire <- missing_parents %in% rel_df$Sire
    is_dam  <- missing_parents %in% rel_df$Dam
    sex_infer <- ifelse(is_sire & !is_dam, 1L, ifelse(is_dam & !is_sire, 2L, 0L))
    rel_df <- bind_rows(
      rel_df,
      tibble(
        ID = missing_parents,
        MatingNumber = NA_character_,
        Sex = NA_character_,
        Dam = NA_real_,
        Sire = NA_real_,
        Sex_k2 = sex_infer
      )
    )
  }
  
  # Prepare vectors for pedigree(); convert 0/"0"/"" to NA
  id_chr <- as.character(rel_df$ID)
  dad_chr <- ifelse(is.na(rel_df$Sire) | rel_df$Sire == 0, NA_character_, as.character(rel_df$Sire))
  mom_chr <- ifelse(is.na(rel_df$Dam)  | rel_df$Dam  == 0, NA_character_, as.character(rel_df$Dam))
  
  # Keep only rows with non-missing IDs
  keep <- !is.na(id_chr) & nzchar(id_chr)
  id_chr  <- id_chr[keep]
  dad_chr <- dad_chr[keep]
  mom_chr <- mom_chr[keep]
  sex_vec <- rel_df$Sex_k2[keep]
  
  # pedigree() requires all parents either NA or present in id list
  ok_parents <- function(par, ids) {
    is.na(par) | par %in% ids
  }
  # If any parent still not present, add them as founders
  extra_parents <- unique(c(dad_chr[!is.na(dad_chr) & !ok_parents(dad_chr, id_chr)],
                            mom_chr[!is.na(mom_chr) & !ok_parents(mom_chr, id_chr)]))
  if (length(extra_parents) > 0) {
    id_chr  <- c(id_chr,  extra_parents)
    dad_chr <- c(dad_chr, rep(NA_character_, length(extra_parents)))
    mom_chr <- c(mom_chr, rep(NA_character_, length(extra_parents)))
    sex_vec <- c(sex_vec, rep(0L, length(extra_parents)))
  }
  
  # Final safety: coerce any leftover "0" to NA
  dad_chr[dad_chr %in% c("0", "")] <- NA_character_
  mom_chr[mom_chr %in% c("0", "")] <- NA_character_
  
  ped_fix <- fixParents(id = id_chr, dadid = dad_chr, momid = mom_chr, sex = sex_vec)
  ped_obj <- pedigree(id = ped_fix$id, dadid = ped_fix$dadid, momid = ped_fix$momid, sex = ped_fix$sex)
  kin <- kinship(ped_obj)
  
  dn <- as.character(dam_id); sn <- as.character(sire_id)
  if (!all(c(dn, sn) %in% rownames(kin))) return(NA_real_)
  round(200 * kin[dn, sn], 1)
}

# ---------------------------- LOAD DATA --------------------------------------
IND_raw <- read.xlsx(IND_PATH, detectDates = TRUE) %>% clean_column_names()
MC_raw  <- read.xlsx(MATING_PATH) %>% clean_column_names()
DNAm    <- read.xlsx(DNAM_PATH) %>% clean_column_names() %>%
  dplyr::select(Species, ExternalSampleName, Age, Sex, EAA)

# -------------------------- COLUMN PICKERS -----------------------------------
IND_cols <- list(
  stock = pick_col(IND_raw, c("STOCK","Stock","Species","strain","Strain","Line","Colony"), required = FALSE, label = "IND.STOCK"),
  id    = pick_col(IND_raw, c("ID","AnimalID","MouseID","IndID","IndividualID"), required = TRUE, label = "IND.ID"),
  sex   = pick_col(IND_raw, c("Sex","SEX","sex"), required = TRUE, label = "IND.Sex"),
  mating= pick_col(IND_raw, c("MatingNumber","Mating","MatingNo","MatingCage","Cage","CageNumber"), required = TRUE, label = "IND.MatingNumber"),
  bday  = pick_col(IND_raw, c("Birthday","DOB","BirthDate"), required = FALSE, label = "IND.Birthday")
)

MC_cols <- list(
  stock = pick_col(MC_raw, c("STOCK","Stock","Species","strain","Strain"), required = FALSE, label = "MC.STOCK"),
  mating= pick_col(MC_raw, c("MatingNumber","Mating","MatingNo","MatingCage","Cage","CageNumber"), required = TRUE, label = "MC.MatingNumber"),
  dam   = pick_col(MC_raw, c("Dam","Mother","DamID"), required = TRUE,  label = "MC.Dam"),
  sire  = pick_col(MC_raw, c("Sire","Father","SireID"), required = TRUE, label = "MC.Sire"),
  date  = pick_col(MC_raw, c("DateofMating","MatingDate","Date"), required = FALSE, label = "MC.DateofMating")
)

# -------------------------- MAIN PIPELINE ------------------------------------


for (GENERATION_DEPTH in c(1:7)){
for (species in TARGET_SPECIES) {
  sp <- toupper(species)
  
  # IND subset
  IND <- IND_raw %>%
    transmute(
      STOCK           = if (!is.na(IND_cols$stock)) .data[[IND_cols$stock]] else NA_character_,
      ID_raw          = .data[[IND_cols$id]],
      Sex             = .data[[IND_cols$sex]],
      Mating_raw      = .data[[IND_cols$mating]],
      Birthday        = if (!is.na(IND_cols$bday)) .data[[IND_cols$bday]] else NA
    ) %>%
    mutate(
      STOCK = normalize_stock(STOCK),
      # Remove species code from ID & MatingNumber BEFORE normalization (critical)
      ID_str          = stringr::str_remove_all(as.character(ID_raw),  regex(sp, ignore_case = TRUE)),
      Mating_str      = stringr::str_remove_all(as.character(Mating_raw), regex(sp, ignore_case = TRUE)),
      ID              = normalize_id(ID_str),
      MatingNumber = clean_mating(Mating_str),
      Sex             = toupper(gsub("[^[:alnum:]]", "", as.character(Sex))),
      Sex = case_when(
        Sex %in% c("F","FEMALE","FEM","FEMALES") ~ "F",
        Sex %in% c("M","MALE","MALES") ~ "M",
        TRUE ~ Sex
      ),
      Birthday = suppressWarnings(as.Date(Birthday, origin = "1899-12-30")),
      BirthMonth = suppressWarnings(lubridate::month(Birthday)),
      BirthYear  = suppressWarnings(lubridate::year(Birthday))
    ) %>%
    filter(is.na(BirthYear) | BirthYear <= 2024) %>%
    filter(!is.na(ID)) %>%
    { if (!is.na(IND_cols$stock)) filter(., str_detect(STOCK, fixed(sp))) else . } %>%
    distinct(ID, .keep_all = TRUE)
  
  # DAM/SIRE subset
  DAMSIRE <- MC_raw %>%
    transmute(
      STOCK        = if (!is.na(MC_cols$stock)) .data[[MC_cols$stock]] else NA_character_,
      Mating_raw   = .data[[MC_cols$mating]],
      Dam_raw      = .data[[MC_cols$dam]],
      Sire_raw     = .data[[MC_cols$sire]],
      DateofMating = if (!is.na(MC_cols$date)) .data[[MC_cols$date]] else NA
    ) %>%
    mutate(
      STOCK        = normalize_stock(STOCK),
      # Remove species code from all three before cleaning/normalizing
      Mating_str   = stringr::str_remove_all(as.character(Mating_raw), regex(sp, ignore_case = TRUE)),
      Dam_str      = stringr::str_remove_all(as.character(Dam_raw),    regex(sp, ignore_case = TRUE)),
      Sire_str     = stringr::str_remove_all(as.character(Sire_raw),   regex(sp, ignore_case = TRUE)),
      MatingNumber = clean_mating(Mating_str),
      Dam          = normalize_id(Dam_str),
      Sire         = normalize_id(Sire_str),
      Dam          = na_if(Dam, 0),
      Sire         = na_if(Sire, 0),
      DateofMating = suppressWarnings(as.Date(DateofMating, origin = "1899-12-30"))
    ) %>%
    { if (!is.na(MC_cols$stock)) filter(., str_detect(STOCK, fixed(sp))) else . } %>%
    filter(!is.na(MatingNumber)) %>%
    distinct(MatingNumber, .keep_all = TRUE)
  
  ped_data <- IND %>%
    left_join(DAMSIRE, by = "MatingNumber") %>%
    dplyr::select(ID, MatingNumber, Sex, Dam, Sire) %>%
    distinct(ID, .keep_all = TRUE)
  
  # DNAm samples
  DNAm_sp <- DNAm %>%
    filter(toupper(Species) == sp) %>%
    mutate(SampleID = normalize_id(ExternalSampleName))
  
  message(glue("[{species}] IND={nrow(IND)} DAMSIRE={nrow(DAMSIRE)} ped_data={nrow(ped_data)} DNAm={nrow(DNAm_sp)}"))
  if (nrow(IND) == 0) {
    message("No IND rows matched this species. Check STOCK values and column detection.")
  }
  
  res  <- vector("list", nrow(DNAm_sp))
  
  for (i in seq_len(nrow(DNAm_sp))) {
    r   <- DNAm_sp[i, ]
    sid <- r$SampleID
    
    IND_row <- IND %>% filter(ID == sid)
    cage <- if (nrow(IND_row) > 0) IND_row$MatingNumber[1] else NA_character_
    
    ds   <- DAMSIRE %>% filter(MatingNumber == cage)
    dam  <- if (nrow(ds) > 0) ds$Dam[1]  else NA_real_
    sire <- if (nrow(ds) > 0) ds$Sire[1] else NA_real_
    
    if (is.na(dam) && is.na(sire)) {
      pd <- ped_data %>% filter(ID == sid)
      if (nrow(pd) > 0) { dam <- pd$Dam[1]; sire <- pd$Sire[1] }
    }
    
    parental_relatedness <- compute_parent_relatedness(dam, sire, IND, DAMSIRE, ped_data, GENERATION_DEPTH)
    
    res[[i]] <- tibble(
      Species = species,
      ExternalSampleName = r$ExternalSampleName,
      Age = r$Age,
      Sex = r$Sex,
      EAA = r$EAA,
      SampleID = sid,
      Match_ID = if (nrow(IND_row) > 0) IND_row$ID[1] else NA_real_,
      MatingNumber = cage,
      Dam = dam,
      Sire = sire,
      ParentalRelatedness = parental_relatedness 
    )
  }
  
  out_df <- bind_rows(res) %>% arrange(ExternalSampleName)
  
  
  
  # ----------------- HEATMAP BY EAA (w/ ANNOTATION) ------------------
  message(glue("[{species}] Generating heatmap by EAA..."))
  
  # 1. Get ordering data
  ordering_data_crr <- out_df %>%
    filter(SampleID %in% target_ids_present) %>%
    distinct(SampleID, .keep_all = TRUE) %>% 
    arrange(EAA)
  
  # 2. Get ordered IDs (for subsetting) and Names (for labels)
  ordered_ids_crr   <- as.character(ordering_data_crr$SampleID)
  ordered_names_crr <- ordering_data_crr$ExternalSampleName
  
  # 3. Re-order the matrix *by SampleID*
  matrix_ordered_crr <- relatedness_matrix[ordered_ids_crr, ordered_ids_crr]
  
  # 4. Apply the ExternalSampleName labels *after* ordering
  rownames(matrix_ordered_crr) <- ordered_names_crr
  colnames(matrix_ordered_crr) <- ordered_names_crr
  
  # 5. Melt for ggplot2
  matrix_melted_crr <- as.data.frame(matrix_ordered_crr) %>%
    tibble::rownames_to_column("ID1") %>%
    tidyr::pivot_longer(-ID1, names_to = "ID2", values_to = "Relatedness") %>%
    mutate(
      ID1 = factor(ID1, levels = ordered_names_crr),
      ID2 = factor(ID2, levels = ordered_names_crr)
    )
  
  # 6. Create MAIN heatmap plot
  g_heatmap_crr <- ggplot(matrix_melted_crr, aes(x = ID1, y = ID2, fill = Relatedness)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(name = "Relatedness (%)", option = "C") +
    theme_minimal(base_size = 8) +
    labs(
      title = glue("{species} - All-Pairs Sample Relatedness"),
      x = NULL, # X-axis label will be on the annotation plot
      y = "Sample"
    ) +
    theme(
      axis.text.x = element_blank(), # X-axis text removed
      axis.ticks.x = element_blank(), # X-axis ticks removed
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5) # Remove bottom margin
    )
  
  # 7. Create ANNOTATION plot
  plot_data_annot_crr <- ordering_data_crr %>%
    mutate(ExternalSampleName = factor(ExternalSampleName, levels = ordered_names_crr))
  
  g_annot_crr <- ggplot(plot_data_annot_crr, aes(x = ExternalSampleName, y = 1, fill = EAA)) +
    geom_tile(color = "grey50", linewidth = 0.1) +
    # Use a different, nice color scale for the annotation
    scale_fill_distiller(palette = "RdYlBu", name = "EAA") +
    scale_x_discrete(expand = c(0, 0)) + # Remove gaps
    scale_y_continuous(expand = c(0, 0)) + # Remove gaps
    theme_minimal(base_size = 8) +
    labs(
      x = "Sample", 
      y = NULL,
      subtitle = glue("Ordered by EAA (n={length(ordered_names_crr)}, depth={GENERATION_DEPTH} gen)")
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5) # Remove top margin
    )
  
  # 8. STITCH plots together
  g_combined_crr <- g_heatmap_crr / g_annot_crr + 
    plot_layout(heights = c(20, 1)) # Heatmap is 20x taller than bar
  
  # 9. Save combined plot
  heatmap_path_crr <- file.path(OUT_DIR, glue("{species} - DNAm relatedness heatmap by EAA (gen {GENERATION_DEPTH}).png"))
  # Adjust height to make room for the annotation bar
  ggsave(heatmap_path_crr, g_combined_crr, width = 10, height = 9, units = "in", dpi = 300)
  message(glue("[{species}] wrote -> {heatmap_path_crr}"))
  
  
  # ------------------- ALL-PAIRS RELATEDNESS MATRIX --------------------
  message(glue("[{species}] Calculating all-pairs relatedness matrix..."))
  
  # Get all unique, non-NA sample IDs from the results
  sample_ids_norm <- stats::na.omit(unique(out_df$SampleID))
  
  # Prepare for matrix & heatmap
  relatedness_matrix <- NULL
  matrix_ordered <- NULL
  
  if (length(sample_ids_norm) < 2) {
    message(glue("[{species}] Skipping matrix/heatmap: < 2 samples with valid IDs."))
  } else {
    
    # 1. Get all ancestors for ALL samples up to the specified depth
    all_ped_ids <- unique(c(
      sample_ids_norm, 
      get_ancestors(sample_ids_norm, ped_data, GENERATION_DEPTH)
    ))
    
    # 2. Build one large pedigree table
    rel_all_df <- IND %>%
      filter(ID %in% all_ped_ids) %>%
      distinct(ID, .keep_all = TRUE) %>%
      left_join(DAMSIRE, by = "MatingNumber") %>%
      distinct(ID, .keep_all = TRUE) %>%
      mutate(
        Dam  = na_if(Dam, 0),
        Sire = na_if(Sire, 0),
        Sex_k2 = case_when(
          toupper(Sex) %in% c("M","MALE") ~ 1L,
          toupper(Sex) %in% c("F","FEMALE") ~ 2L,
          TRUE ~ 0L
        )
      )
    
    # Ensure all focal samples exist
    missing_samples <- setdiff(sample_ids_norm, rel_all_df$ID)
    if (length(missing_samples) > 0) {
      sample_sex_info <- out_df %>% 
        filter(SampleID %in% missing_samples) %>%
        distinct(SampleID, Sex) %>%
        mutate(Sex_k2 = case_when(
          toupper(Sex) %in% c("M","MALE") ~ 1L,
          toupper(Sex) %in% c("F","FEMALE") ~ 2L,
          TRUE ~ 0L
        ))
      
      rel_all_df <- bind_rows(
        rel_all_df,
        tibble(
          ID = sample_sex_info$SampleID,
          Sex = sample_sex_info$Sex,
          Sex_k2 = sample_sex_info$Sex_k2,
          MatingNumber = NA_character_, Dam = NA_real_, Sire = NA_real_
        )
      )
    }
    
    # Ensure *all* referenced parents exist as rows
    all_parent_ids <- unique(stats::na.omit(c(rel_all_df$Dam, rel_all_df$Sire)))
    missing_parents <- setdiff(all_parent_ids, rel_all_df$ID)
    if (length(missing_parents) > 0) {
      is_sire <- missing_parents %in% rel_all_df$Sire
      is_dam  <- missing_parents %in% rel_all_df$Dam 
      sex_infer <- ifelse(is_sire & !is_dam, 1L, ifelse(is_dam & !is_sire, 2L, 0L))
      rel_all_df <- bind_rows(
        rel_all_df,
        tibble(
          ID = missing_parents,
          MatingNumber = NA_character_, Sex = NA_character_,
          Dam = NA_real_, Sire = NA_real_, Sex_k2 = sex_infer
        )
      )
    }
    
    # 3. Prepare vectors for pedigree()
    id_chr  <- as.character(rel_all_df$ID)
    dad_chr <- ifelse(is.na(rel_all_df$Sire) | rel_all_df$Sire == 0, NA_character_, as.character(rel_all_df$Sire))
    mom_chr <- ifelse(is.na(rel_all_df$Dam)  | rel_all_df$Dam  == 0, NA_character_, as.character(rel_all_df$Dam))
    
    keep <- !is.na(id_chr) & nzchar(id_chr)
    id_chr  <- id_chr[keep]
    dad_chr <- dad_chr[keep]
    mom_chr <- mom_chr[keep]
    sex_vec <- rel_all_df$Sex_k2[keep]
    
    ok_parents <- function(par, ids) { is.na(par) | par %in% ids }
    extra_parents <- unique(c(
      dad_chr[!is.na(dad_chr) & !ok_parents(dad_chr, id_chr)],
      mom_chr[!is.na(mom_chr) & !ok_parents(mom_chr, id_chr)]
    ))
    if (length(extra_parents) > 0) {
      id_chr  <- c(id_chr,  extra_parents)
      dad_chr <- c(dad_chr, rep(NA_character_, length(extra_parents)))
      mom_chr <- c(mom_chr, rep(NA_character_, length(extra_parents)))
      sex_vec <- c(sex_vec, rep(0L, length(extra_parents)))
    }
    
    dad_chr[dad_chr %in% c("0", "")] <- NA_character_
    mom_chr[mom_chr %in% c("0", "")] <- NA_character_
    
    # 4. Create pedigree and kinship matrix
    ped_fix <- fixParents(id = id_chr, dadid = dad_chr, momid = mom_chr, sex = sex_vec)
    ped_obj_all <- pedigree(id = ped_fix$id, dadid = ped_fix$dadid, momid = ped_fix$momid, sex = ped_fix$sex)
    kin_all <- kinship(ped_obj_all)
    
    # 5. Subset matrix to *only* our samples
    target_ids_chr <- as.character(sample_ids_norm)
    target_ids_present <- target_ids_chr[target_ids_chr %in% rownames(kin_all)]
    
    if (length(target_ids_present) < 2) {
      message(glue("[{species}] Skipping matrix/heatmap: < 2 samples found in final pedigree."))
    } else {
      # Keep the matrix indexed by SampleID (which is unique)
      kin_subset <- kin_all[target_ids_present, target_ids_present]
      relatedness_matrix <- round(200 * kin_subset, 1)
      
      # ----------------- HEATMAP BY PARENTAL RELATEDNESS (w/ ANNOTATION) -----------------
      message(glue("[{species}] Generating heatmap by parental relatedness..."))
      
      # 1. Get ordering data
      ordering_data_pr <- out_df %>%
        filter(SampleID %in% target_ids_present) %>%
        distinct(SampleID, .keep_all = TRUE) %>%
        arrange(is.na(ParentalRelatedness), ParentalRelatedness)
      
      # 2. Get ordered IDs and Names
      ordered_ids_pr    <- as.character(ordering_data_pr$SampleID)
      ordered_names_pr <- ordering_data_pr$ExternalSampleName
      
      # 3. Re-order the matrix *by SampleID*
      matrix_ordered_pr <- relatedness_matrix[ordered_ids_pr, ordered_ids_pr]
      
      # 4. Apply labels *after* ordering
      rownames(matrix_ordered_pr) <- ordered_names_pr
      colnames(matrix_ordered_pr) <- ordered_names_pr
      
      # 5. Melt
      matrix_melted_pr <- as.data.frame(matrix_ordered_pr) %>%
        tibble::rownames_to_column("ID1") %>%
        tidyr::pivot_longer(-ID1, names_to = "ID2", values_to = "Relatedness") %>%
        mutate(
          ID1 = factor(ID1, levels = ordered_names_pr),
          ID2 = factor(ID2, levels = ordered_names_pr)
        )
      
      # 6. Create MAIN heatmap plot
      g_heatmap_pr <- ggplot(matrix_melted_pr, aes(x = ID1, y = ID2, fill = Relatedness)) +
        geom_tile(color = "white", linewidth = 0.1) +
        scale_fill_viridis_c(name = "Relatedness (%)", option = "C") +
        theme_minimal(base_size = 8) +
        labs(
          title = glue("{species} - All-Pairs Sample Relatedness"),
          x = NULL, # X-axis label will be on the annotation plot
          y = "Sample"
        ) +
        theme(
          axis.text.x = element_blank(), # X-axis text removed
          axis.ticks.x = element_blank(), # X-axis ticks removed
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5) # Remove bottom margin
        )
      
      # 7. Create ANNOTATION plot
      plot_data_annot_pr <- ordering_data_pr %>%
        mutate(ExternalSampleName = factor(ExternalSampleName, levels = ordered_names_pr))
      
      g_annot_pr <- ggplot(plot_data_annot_pr, aes(x = ExternalSampleName, y = 1, fill = ParentalRelatedness)) +
        geom_tile(color = "grey50", linewidth = 0.1) +
        # Use a different, nice color scale for the annotation
        scale_fill_distiller(palette = "PuOr", name = "Parental\nRelatedness") +
        scale_x_discrete(expand = c(0, 0)) + # Remove gaps
        scale_y_continuous(expand = c(0, 0)) + # Remove gaps
        theme_minimal(base_size = 8) +
        labs(
          x = "Sample", 
          y = NULL,
          subtitle = glue("Ordered by Parental Relatedness (n={length(ordered_names_pr)}, depth={GENERATION_DEPTH} gen)")
        ) +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5) # Remove top margin
        )
      
      # 8. STITCH plots together
      g_combined_pr <- g_heatmap_pr / g_annot_pr + 
        plot_layout(heights = c(20, 1)) # Heatmap is 20x taller than bar
      
      # 9. Save combined plot
      heatmap_path_pr <- file.path(OUT_DIR, glue("{species} - DNAm relatedness heatmap by ParentalRelatedness (gen {GENERATION_DEPTH}).png"))
      # Adjust height to make room for the annotation bar
      ggsave(heatmap_path_pr, g_combined_pr, width = 10, height = 9, units = "in", dpi = 300)
      message(glue("[{species}] wrote -> {heatmap_path_pr}"))
      
      # Set the main matrix_ordered variable for saving to Excel (USING PR ORDER)
      matrix_ordered <- matrix_ordered_pr
    }
  }
  
  # ------------------------- UPDATED EXCEL OUTPUT ------------------------
  
  out_path <- file.path(OUT_DIR, glue(OUT_FILE_TMPL, species = species, gen_depth = GENERATION_DEPTH))
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add the original results
  addWorksheet(wb, "SampleParentalRelatedness")
  writeData(wb, "SampleParentalRelatedness", out_df)
  
  # Add the new matrix if it was created
  if (!is.null(matrix_ordered)) {
    addWorksheet(wb, "AllPairsMatrix_by_PR")
    # This saves the matrix ordered by Parental Relatedness
    matrix_df <- as.data.frame(matrix_ordered) %>%
      tibble::rownames_to_column("Sample")
    writeData(wb, "AllPairsMatrix_by_PR", matrix_df, rowNames = FALSE)
  }
  
  # Save the workbook
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message(glue("[{species}] wrote -> {out_path}"))
  
} # End of the main "for (species ...)" loop
}

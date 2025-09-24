

suppressPackageStartupMessages({
  library(openxlsx)
  library(tidyverse)
  library(lubridate)
  library(glue)
  library(kinship2)
})

# ---------------------------- CONFIG -----------------------------------------
IND_PATH        <- "./data/Peromyscus.xlsx"
MATING_PATH     <- "./data/Mating Records.xlsx"
DNAM_PATH       <- "./data/DNAm mice info -BW dataset.xlsx"
TARGET_SPECIES  <- c("BW")
MIN_GENERATIONS <- 5
MAX_GENERATIONS <- 10
OUT_DIR         <- "./relatedness"
OUT_FILE_TMPL   <- "{species} - DNAm parental relatedness (gens {min_g}-{max_g}).xlsx"

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

compute_parent_relatedness <- function(dam_id, sire_id, IND, DAMSIRE, ped_data, depth) {
  if (any(is.na(c(dam_id, sire_id)))) return(NA_real_)
  
  # Scope by ancestors up to depth
  base_ids <- unique(c(dam_id, sire_id))
  rel_ids  <- unique(c(base_ids, get_ancestors(base_ids, ped_data, depth)))
  
  # Start from IND rows for those IDs
  rel <- IND %>%
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
  if (!(dam_id %in% rel$ID)) rel <- bind_rows(rel, tibble(ID = dam_id, MatingNumber = NA_character_, Sex = "F", Dam = NA_real_, Sire = NA_real_, Sex_k2 = 2L))
  if (!(sire_id %in% rel$ID)) rel <- bind_rows(rel, tibble(ID = sire_id, MatingNumber = NA_character_, Sex = "M", Dam = NA_real_, Sire = NA_real_, Sex_k2 = 1L))
  
  # Ensure *all* referenced parents exist as rows (founders if absent)
  all_parent_ids <- unique(stats::na.omit(c(rel$Dam, rel$Sire)))
  missing_parents <- setdiff(all_parent_ids, rel$ID)
  if (length(missing_parents) > 0) {
    # Infer sex where possible from being listed as Sire/Dam
    is_sire <- missing_parents %in% rel$Sire
    is_dam  <- missing_parents %in% rel$Dam
    sex_infer <- ifelse(is_sire & !is_dam, 1L, ifelse(is_dam & !is_sire, 2L, 0L))
    rel <- bind_rows(
      rel,
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
  id_chr <- as.character(rel$ID)
  dad_chr <- ifelse(is.na(rel$Sire) | rel$Sire == 0, NA_character_, as.character(rel$Sire))
  mom_chr <- ifelse(is.na(rel$Dam)  | rel$Dam  == 0, NA_character_, as.character(rel$Dam))
  
  # Keep only rows with non-missing IDs
  keep <- !is.na(id_chr) & nzchar(id_chr)
  id_chr  <- id_chr[keep]
  dad_chr <- dad_chr[keep]
  mom_chr <- mom_chr[keep]
  sex_vec <- rel$Sex_k2[keep]
  
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
  dplyr::select(Species, ExternalSampleName, Age, Sex)

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
for (species in TARGET_SPECIES) {
  sp <- toupper(species)
  
  # IND subset
  IND <- IND_raw %>%
    transmute(
      STOCK        = if (!is.na(IND_cols$stock)) .data[[IND_cols$stock]] else NA_character_,
      ID_raw       = .data[[IND_cols$id]],
      Sex          = .data[[IND_cols$sex]],
      Mating_raw   = .data[[IND_cols$mating]],
      Birthday     = if (!is.na(IND_cols$bday)) .data[[IND_cols$bday]] else NA
    ) %>%
    mutate(
      STOCK = normalize_stock(STOCK),
      # Remove species code from ID & MatingNumber BEFORE normalization (critical)
      ID_str       = stringr::str_remove_all(as.character(ID_raw),  regex(sp, ignore_case = TRUE)),
      Mating_str   = stringr::str_remove_all(as.character(Mating_raw), regex(sp, ignore_case = TRUE)),
      ID           = normalize_id(ID_str),
      MatingNumber = clean_mating(Mating_str),
      Sex          = toupper(gsub("[^[:alnum:]]", "", as.character(Sex))),
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
  
  GENS <- seq.int(MIN_GENERATIONS, MAX_GENERATIONS)
  res  <- vector("list", nrow(DNAm_sp))
  
  for (i in seq_len(nrow(DNAm_sp))) {
    r   <- DNAm_sp[i, ]
    sid <- r$SampleID
    
    IND_row <- IND %>% filter(ID == sid)
    cage <- if (nrow(IND_row) > 0) IND_row$MatingNumber[1] else NA_character_
    
    ds    <- DAMSIRE %>% filter(MatingNumber == cage)
    dam   <- if (nrow(ds) > 0) ds$Dam[1]  else NA_real_
    sire  <- if (nrow(ds) > 0) ds$Sire[1] else NA_real_
    
    if (is.na(dam) && is.na(sire)) {
      pd <- ped_data %>% filter(ID == sid)
      if (nrow(pd) > 0) { dam <- pd$Dam[1]; sire <- pd$Sire[1] }
    }
    
    rel_cols <- map_dbl(GENS, ~ compute_parent_relatedness(dam, sire, IND, DAMSIRE, ped_data, .x))
    
    res[[i]] <- tibble(
      Species = species,
      ExternalSampleName = r$ExternalSampleName,
      Age = r$Age,
      Sex = r$Sex,
      SampleID = sid,
      Match_ID = if (nrow(IND_row) > 0) IND_row$ID[1] else NA_real_,
      MatingNumber = cage,
      Dam = dam,
      Sire = sire,
      !!!setNames(as.list(rel_cols), paste0("Relatedness_g", GENS))
    )
  }
  
  out_df <- bind_rows(res) %>% arrange(ExternalSampleName)
  
  out_path <- file.path(OUT_DIR, glue(OUT_FILE_TMPL, species = species, min_g = MIN_GENERATIONS, max_g = MAX_GENERATIONS))
  write.xlsx(out_df, file = out_path, overwrite = TRUE)
  
  message(glue("[{species}] wrote -> {out_path}"))
}

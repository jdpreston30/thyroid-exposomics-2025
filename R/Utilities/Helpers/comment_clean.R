#! Helper
#' Clean Hierarchical Header Comments Across All Pipeline Scripts
#'
#' For every `.R` file under \code{root_dir}, applies three normalizations to
#' lines that begin with \code{#*}, \code{#+}, \code{#-}, \code{#_}, or
#' \code{#!}:
#'
#' 1. **Renumber**: fixes the leading numeric segment of every \code{#+} and
#'    \code{#-} header so that (a) the top-level number matches the file prefix
#'    (e.g. \code{06_*.R} → prefix \code{6}; \code{00e_*.R} → prefix \code{0e})
#'    and (b) subsections are always contiguous starting at \code{1} — gaps are
#'    closed and \code{.0} starts are shifted up (e.g. if a \code{.1} subsection
#'    is deleted, the remaining \code{.2}, \code{.3}, ... become \code{.1},
#'    \code{.2}, ...). Only \code{#*}, \code{#+}, \code{#-} tagged lines are
#'    touched; \code{#!} and \code{#_} are left unchanged.
#' 2. **Format**: enforces exactly \code{#tag<space>number:<space>message} —
#'    inserts a missing colon, removes extra spaces around the colon, normalizes
#'    to exactly one space between tag and number.
#' 3. **Tag-aware casing**: `#*` and `#+` messages are Title-Cased (each word
#'    capitalized; small words like "and"/"of"/"the" kept lowercase unless they
#'    are the first word), while `#-` messages are sentence-cased (only the first
#'    word capitalized). Both preserve all-caps acronyms (IARC, GHS), mixed-case
#'    names (PubChem, ClassyFire), tokens with digits/punctuation (ST3,
#'    o-Toluidine, ~3%), and the protected proper-noun list. `#!` and `#_` are
#'    never re-cased.
#'
#' Only header lines are touched — all code lines pass through unchanged.
#'
#' @param root_dir Character. Root directory to scan recursively. Default
#'   \code{"R/Scripts"}.
#' @param dry_run Logical. If \code{TRUE}, prints diffs without writing files.
#'   Default \code{FALSE}.
#' @param verbose Logical. If \code{TRUE}, prints per-file change counts.
#'   Default \code{TRUE}.
#'
#' @return Invisibly returns a named list: \code{files_scanned},
#'   \code{files_changed}, \code{line_changes}.
#'
#' @export
comment_clean <- function(root_dir = "R/Scripts", dry_run = FALSE, verbose = TRUE) {
  if (!dir.exists(root_dir)) stop("comment_clean(): root_dir does not exist: ", root_dir)

  r_files <- list.files(root_dir, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  if (length(r_files) == 0) {
    if (isTRUE(verbose)) message("comment_clean(): no .R files found under ", root_dir)
    return(invisible(list(files_scanned = 0L, files_changed = 0L, line_changes = 0L)))
  }

  # --- inner helpers --------------------------------------------------------

  extract_file_prefix <- function(path) {
    reg <- regmatches(basename(path), regexec("^([0-9]+)([A-Za-z]?)", basename(path)))[[1]]
    if (length(reg) < 3) return(NA_character_)
    paste0(as.integer(reg[2]), reg[3])   # "00" -> "0", "00e" -> "0e", "11" -> "11"
  }

  # 1. Sequential renumber: fix file-level prefix AND make #+/#- subsections
  # contiguous starting at 1 with no gaps. Must operate on all lines at once
  # (single-pass sequential counter), not line-by-line.
  renumber_subsections <- function(lines, file_prefix) {
    pat <- "^(\\s*#([*+\\-])\\s+)[0-9]+[A-Za-z]?(?:\\.[0-9]+)*(.*)$"
    plus_counter  <- 0L
    minus_counter <- 0L
    result <- lines
    for (i in seq_along(lines)) {
      m   <- regexec(pat, lines[i])
      reg <- regmatches(lines[i], m)[[1]]
      if (length(reg) == 0) next
      tag_prefix <- reg[2]
      tag_char   <- reg[3]
      rest       <- sub("^[A-Za-z](?=:)", "", reg[4], perl = TRUE)   # strip letter suffixes (14.8a → 14.8)
      if (tag_char == "*") {
        plus_counter  <- 0L
        minus_counter <- 0L
        result[i] <- paste0(tag_prefix, file_prefix, rest)
      } else if (tag_char == "+") {
        plus_counter  <- plus_counter + 1L
        minus_counter <- 0L
        result[i] <- paste0(tag_prefix, file_prefix, ".", plus_counter, rest)
      } else if (tag_char == "-") {
        if (plus_counter == 0L) next
        minus_counter <- minus_counter + 1L
        result[i] <- paste0(tag_prefix, file_prefix, ".", plus_counter, ".", minus_counter, rest)
      }
    }
    result
  }

  # 2. Tag-aware casing:
  #    - #* / #+ headers -> Title Case (each word capitalized; small words lower)
  #    - #-  headers     -> sentence case (only the first word capitalized)
  #    Applied only to numbered headers; #! and #_ are never re-cased.
  # Words in .proper_nouns are always preserved regardless of capitalization.
  .proper_nouns <- c(
    # software / names / stats
    "Arial", "PostScript", "GitHub", "RStudio", "PubChem",
    "Excel", "Docker", "Emory", "Tukey", "Fisher", "Kruskal", "Wallis", "Restek",
    # real directory names (keep capitalized even mid-message)
    "Components", "Supplementary", "Sections", "References",
    # DTC type names
    "Papillary", "Follicular",
    # single-word chemical names (multi-word / all-caps / digit forms are already
    # preserved by .is_special; these plain Title-Case names are not)
    "Methylparaben", "Cyfluthrin", "Napropamide", "Molinate", "Atrazine",
    "Bupirimate", "Resmethrin", "Vernolate", "Chrysene", "Naphthalene",
    "Anthracene", "Fluoranthene", "Pyrene", "Acetophenone", "Ethylan",
    "Metalaxyl", "Terbuthylazine", "Prosulfuron", "Flucythrinate", "Menthone",
    "Pentachlorophenol", "Phenacetin", "Hexanoic", "Dibutyl", "Ethyl"
  )
  # Title-Case small words kept lowercase unless they are the first word.
  .small_words <- c(
    "a", "an", "and", "or", "of", "the", "to", "in", "for", "with",
    "vs", "v", "by", "from", "on", "at", "as", "but", "nor", "per", "via"
  )
  # A word is preserved verbatim (never re-cased) if it is a protected proper noun,
  # contains a digit, is an ALL-CAPS acronym (>=2 letters), has an internal
  # uppercase letter (mixed case), or contains punctuation other than a hyphen.
  .is_special <- function(w) {
    w %in% .proper_nouns ||
      grepl("[0-9]", w) ||
      grepl("^[A-Z]{2,}$", w) ||
      grepl("^.+[A-Z]", w) ||
      grepl("^[A-Za-z]-", w) ||   # single-letter locant/prefix: o-, m-, N-, z-, t- …
      grepl("[^A-Za-z-]", w)
  }
  .cap_first <- function(w) {
    if (nchar(w) == 0) return(w)
    paste0(toupper(substr(w, 1, 1)), substr(w, 2, nchar(w)))
  }
  # Sentence case (#-): first word capitalized, others lowercased unless special.
  to_sentence_case <- function(text) {
    words <- strsplit(text, " ", fixed = TRUE)[[1]]
    if (length(words) == 0) return(text)
    words[-1] <- vapply(words[-1], function(w) {
      if (.is_special(w)) w else tolower(w)
    }, character(1))
    # Capitalize the first word ONLY if it is not a special token — this keeps
    # lowercase locants (o-, N-) and Greek letters (γ-) intact, e.g. o-Toluidine,
    # γ-BHC are left as-is rather than mangled to O-Toluidine / Γ-BHC.
    if (!.is_special(words[1])) words[1] <- .cap_first(words[1])
    paste(words, collapse = " ")
  }
  # Title Case (#* / #+): capitalize each word; small words stay lower (except
  # the first); special tokens preserved verbatim.
  to_title_case <- function(text) {
    words <- strsplit(text, " ", fixed = TRUE)[[1]]
    if (length(words) == 0) return(text)
    for (i in seq_along(words)) {
      w <- words[i]
      if (.is_special(w)) next
      if (i > 1 && tolower(w) %in% .small_words) {
        words[i] <- tolower(w)
      } else {
        words[i] <- .cap_first(w)
      }
    }
    paste(words, collapse = " ")
  }

  # 3. Normalize format: one space after tag, colon + one space before message
  normalize_format <- function(line) {
    m <- regexec("^(\\s*#([-*+_!]))(\\s+)(\\S.*)?$", line)  # '-' first: literal, not a range
    reg <- regmatches(line, m)[[1]]
    if (length(reg) == 0) return(line)
    tag      <- reg[2]
    tag_char <- reg[3]
    rest     <- if (length(reg) >= 5 && !is.na(reg[5])) reg[5] else ""

    if (tag_char %in% c("*", "+", "-")) {
      # Numbered header: parse number token then message
      m2   <- regexec("^([0-9][^:[:blank:]]*)[[:blank:]]*:?[[:blank:]]*(.*?)[[:blank:]]*$", rest)
      reg2 <- regmatches(rest, m2)[[1]]
      if (length(reg2) >= 3) {
        msg <- if (tag_char %in% c("*", "+")) to_title_case(reg2[3]) else to_sentence_case(reg2[3])
        return(paste0(tag, " ", reg2[2], ": ", msg))
      }
    }

    # #! and #_ : normalize spacing only — no sentence case (may contain code/emphasis)
    paste0(tag, " ", trimws(rest))
  }

  # --- main loop ------------------------------------------------------------

  files_changed <- 0L
  line_changes  <- 0L

  for (f in r_files) {
    file_prefix <- extract_file_prefix(f)
    if (is.na(file_prefix)) next

    lines <- readLines(f, warn = FALSE)
    if (length(lines) == 0) next

    lines_renumbered <- renumber_subsections(lines, file_prefix)
    new_lines <- vapply(lines_renumbered, normalize_format, character(1))

    changed_idx <- which(lines != new_lines)

    if (length(changed_idx) > 0) {
      files_changed <- files_changed + 1L
      line_changes  <- line_changes + length(changed_idx)

      if (isTRUE(verbose)) {
        message(sprintf("comment_clean(): %s | %d line(s) %s",
                        f, length(changed_idx),
                        if (isTRUE(dry_run)) "would change" else "changed"))
        if (isTRUE(dry_run)) {
          for (i in changed_idx) {
            message(sprintf("  line %d  before: %s", i, lines[i]))
            message(sprintf("  line %d  after:  %s", i, new_lines[i]))
          }
        }
      }

      if (!isTRUE(dry_run)) writeLines(new_lines, con = f, useBytes = TRUE)
    }
  }

  if (isTRUE(verbose)) {
    message(sprintf("comment_clean(): scanned=%d, files_%s=%d, line_changes=%d",
                    length(r_files),
                    if (isTRUE(dry_run)) "to_change" else "changed",
                    files_changed, line_changes))
  }

  invisible(list(files_scanned = length(r_files),
                 files_changed = files_changed,
                 line_changes  = line_changes))
}

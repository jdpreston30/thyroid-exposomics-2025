#' Normalize chemical nomenclature for display
#'
#' Applies IUPAC/IARC display conventions to chemical names read from the library
#' and metadata sheets. Corrections live here rather than in the source workbooks
#' because those sit on OneDrive and are not editable from the pipeline.
#'
#' @param x Character vector of chemical names.
#'
#' @return Character vector, same length, with display corrections applied.
#'
fix_chem_nomenclature <- function(x) {
#! Fusion locants take square brackets (IUPAC/IARC: benzo[a]pyrene). Guarded on a letter before "(" and a lowercase letter after ")" so only fusion positions match -- substituent parentheses such as bis(2,3-dibromopropyl), 4,4'-Methylenebis(2-chloroaniline), 3(2H)-furanone and "(linear)" all fail one of the two guards and are left alone.
  x <- stringr::str_replace_all(x, "(?<=[A-Za-z])\\(([a-z](?:,[a-z0-9]+)*)\\)(?=[a-z])", "[\\1]")
#! Cresyl phosphate isomers: o/m/p are ortho/meta/para and must be lowercase. Capital O- and P- are themselves real locants (oxygen, phosphorus), so the source casing names the wrong substitution. CAS confirms: 78-30-8 ortho, 563-04-2 meta, 78-32-0 para.
  x <- stringr::str_replace(x, "^Tri-([OMP])-cresyl", \(m) paste0("Tri-", tolower(stringr::str_sub(m, 5, 5)), "-cresyl"))
#! Structural prefix, always lowercase (CP2351, CAS 140-66-9 = 4-tert-octylphenol).
  x <- stringr::str_replace(x, "-Tert-", "-tert-")
#! "Fluro" is not a chemical prefix; the element is fluorine, so the prefix is "fluoro". The trailing hyphen in "phenoxy-benzoic" is also non-standard (CP2329, CAS 77279-89-1 = 4-fluoro-3-phenoxybenzoic acid).
  x <- stringr::str_replace(x, "Fluro-3-phenoxy-benzoic", "Fluoro-3-phenoxybenzoic")
#! 2-ethylhexyl is one substituent, so its locant must be parenthesised or the "2-" reads as a phthalate
#! ring position. CAS 4376-20-9 = mono(2-ethylhexyl) phthalate; the sibling entries bis(2-Ethylhexyl)phthalate
#! and Tris(2-ethylhexyl) phosphate already parenthesise the same group.
  x <- stringr::str_replace(x, "Mono-2-ethylhexyl phthalate", "Mono(2-ethylhexyl) phthalate")
#! Primed locants take U+2032 PRIME, not an apostrophe (4,4'-Diaminodiphenylmethane). Escaped rather than written literally: a literal prime is read as U+FFFD when the session locale is not UTF-8, which silently corrupts the label. Guarded on digit-comma-digit so ordinary apostrophes are untouched.
  x <- stringr::str_replace_all(x, "(?<=\\d),(\\d+)'", ",\\1\\u2032")
  x
}

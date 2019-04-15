# package names from https://bioconductor.org/packages/3.5/data/annotation/

org_pkg <- fread("PASTE TEXT HERE", header = FALSE)
org_pkg <- org_pkg$V1[org_pkg$V2 == 'Bioconductor Package Maintainer']

# get packages
## try http:// if https:// URLs are not supported
BiocManager::install(org_pkg)

package_organisms <- c()

for (package_name in org_pkg) {

    pkg <- crossmeta:::get_biocpack(package_name)
    package_organisms <- c(package_organisms, AnnotationDbi::species(pkg))
}

names(org_pkg) <- package_organisms

devtools::use_data(gpl_bioc, homologene, sources, token, org_pkg, internal = TRUE, overwrite = TRUE)



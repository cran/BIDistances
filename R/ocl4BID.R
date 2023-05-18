.ocl_system_file <-
    function(...)
{
    file <- system.file(package = "BIDistances", ...)
    readLines(file)
}

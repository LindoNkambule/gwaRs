# Copyright (C) Lindokuhle Nkambule
#
# This file is part of gwaRs
#
# gwaRs is free software: you can redistribute it and/or modify it
# under the terms of the MIT as published by
# the Free Software Foundation.
#
# gwaRs is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# MIT License for more details.
#
#  A copy of the MIT is available at
#  http://www.r-project.org/Licenses/
#


#' @import utils
NULL

.onAttach <- function(lib, pkg,...){
  packageStartupMessage(gwaRsWelcomeMessage())
}


gwaRsWelcomeMessage <- function(){

  cat("Thank you for using gwaRs (version ",utils::packageDescription("gwaRs")$Version, ")\n",
      "To access the overall documentation, go to: https://cran.r-project.org/web/packages/gwaRs/vignettes/gwaRs.html\n",
      "Suggestions and bug-reports can be submitted at: https://github.com/LindoNkambule/gwaRs/issues\n",
      "\n",
      "\tTo suppress this message use:\n",
      "\tsuppressPackageStartupMessages(library(gwaRs))\n",
      sep="")
}

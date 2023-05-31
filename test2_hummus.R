devtools::load_all("../hummus_package")

biboo <- define_grn(
  hummus,
  multilayer_f = "a",
  njobs = 5
  )
library(EthSEQ)

ethseq.Analysis(
  target.vcf = "./1000GP_Genotypes100s.vcf",
  out.dir = "./100s",
  model.gds = "./ReferenceModel.gds",
  cores=1,
  verbose=TRUE,
  composite.model.call.rate = 0.99,
  space="3D")

ethseq.Analysis(
  target.vcf = "./1000GP_Genotypes16s.vcf",
  out.dir = "./16s",
  model.gds = "./ReferenceModel.gds",
  cores=1,
  verbose=TRUE,
  composite.model.call.rate = 0.99,
  space="3D")

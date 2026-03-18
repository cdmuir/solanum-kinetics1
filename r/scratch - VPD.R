# summarize VPD for paper, rename, delete
source("r/header.R")

read_rds("processed-data/rh_curves.rds") |> colnames() |> sort()

li6800_svp = function(T_degreeC) {
  0.61365 * exp(17.502 * T_degreeC / (240.97 + T_degreeC))
}

SVPleaf
VPDleaf = (SVPleaf - H2O_s * (Pa + deltaPcham) / 1000)
tmp = read_licor("raw-data/licor/2022-03-14-0923_logdata")
tmp$Pa
tmp$VPDleaf
tmp$Tleaf
tmp$Tair
li6800_svp(tmp$Tleaf)
tmp$SVPleaf
tmp$SVPair


rh_curves = read_rds("data/rh_curves.rds") |> 
  mutate(
    SVPleaf = li6800_svp(Tleaf))
rh_curves |>
  filter(Tleaf < 24) |>
  select(acc_id) |>
  print(n = Inf)

rh_curves |>
  filter(Tleaf > 25.5) |>
  select(acc_id) |>
  print(n = Inf)

rh_curves |>
  colnames() |> sort()
ggplot(filter(rh_curves, acc_id != "LA2744-G"), aes(x = Tleaf)) +
  geom_histogram() +
  facet_wrap(~light_intensity)

rh_curves$Pa

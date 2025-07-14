source(here::here("Summary metrics/activity_bouts.R"))
source(here::here("Summary metrics/activity_index.R"))
source(here::here("Summary metrics/cum_power.R"))
source(here::here("Summary metrics/submovements.R"))

source(here::here("Summarise activity bouts.R"))

# system("dx download project-GgQB2KQJkFy6345bP90G5Q9Z:file-J18jF78JkFy05qK3kG4yGp8B -o Dx/")

set.seed(123)

# Load list of eids into R
Data = readr::read_csv("Dx/Matched_cases_controls_activity.csv")

sample_ids = (Data |>
                group_by(subclass) |>
                filter(n() == 6 & group == "Case") |>
                slice(1))$subclass |>
  sample(10)

Sample_eids = (Data |>
                 filter(subclass %in% sample_ids) |>
                 select(eid))$eid

#Sample_eids = Sample_eids[1:2]

files_raw = paste0("/Bulk/Activity/Raw/", stringr::str_extract(Sample_eids, "^[0-9]{2}"),
                   "/", Sample_eids, "_90001_0_0.cwa")

# Download files associated with sample
# sapply(files_raw, function(x) system(paste0("dx download ", x, " -o CWAs")))

done_files_str = "[0-9][0-9][0-9][0-9][0-9][0-9][0-9]"

done_eids = stringr::str_match(list.files("AI bouts"), done_files_str) |> as.numeric()

# Filtering out files I have already done
Sample_eids = Sample_eids[!(Sample_eids %in% done_eids)]

lapply(Sample_eids, summarise_avtivity_bouts)

# Upload AI files to UKB
# system(dx upload "AI bouts/" -r)

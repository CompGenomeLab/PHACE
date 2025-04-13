using PlmDCA

id = ARGS[1]

output_dir = "/cta/users/nkuru/DCA/Results"
fasta_path = joinpath("/cta/groups/adebali/phact_data_nurdan/MaskedMSA_RAxMLProteins", "$(id)_MaskedMSA.fasta")
txt_file_path = joinpath(output_dir, "$(id)_dca_scores.txt")

plm_out = plmdca(fasta_path)
scores = plm_out.score

open(txt_file_path, "w") do file
    for (i, j, val) in scores
        println(file, "$i $j $val")
    end
end


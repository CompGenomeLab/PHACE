using GaussDCA

id = ARGS[1]

output_dir = "/cta/users/nkuru/GaussDCA/Results"
fasta_path = joinpath("MaskedMSA_RAxMLProteins", "$(id)_MaskedMSA.fasta")
txt_file_path = joinpath(output_dir, "$(id)_dca_scores.txt")

gdca_out = gdca(fasta_path)

# Check!
scores = gdca_out.score

open(txt_file_path, "w") do file
    for (i, j, val) in scores
        println(file, "$i $j $val")
    end
end

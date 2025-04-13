using MIToS.MSA
using MIToS.Information

# Assuming id is passed as a command-line argument
id = ARGS[1]

output_dir = "MutualInfo/Results"
output_dir2 = "MutualInfo/Results2"
fasta_path = joinpath("MaskedMSA_RAxMLProteins", "$(id)_MaskedMSA.fasta")
txt_file_path = joinpath(output_dir, "$(id)_MIp_scores.txt")
txt_file_path2 = joinpath(output_dir2, "$(id)_MIp_scores.txt")

msa = read(fasta_path, FASTA)
ZMIp, MIp = buslje09(msa)
scores = MIp

open(txt_file_path, "w") do file
    for i in 1:size(scores, 1)
        for j in 1:size(scores, 2)
            val = scores[i, j]
            println(file, "$i $j $val")
        end
    end
end



scores = ZMIp

open(txt_file_path2, "w") do file
    for i in 1:size(scores, 1)
        for j in 1:size(scores, 2)
            val = scores[i, j]
            println(file, "$i $j $val")
        end
    end
end


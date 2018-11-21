    for i in `seq 1 22`
    do
    awk '{print $3}' /spin1/users/zhangh24/BCAC/impute_plink_onco/chr$i\_ld_clump.clumped >> /spin1/users/zhangh24/BCAC/impute_plink_onco/clump_snp
    done

tail -n +2 clump_snp > clump_snp_new
rm clump_snp
mv clump_snp_new clump_snp
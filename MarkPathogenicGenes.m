function [newdata,risky_length, risky_met_genes_length,risky_met_genes,before_exp]= MarkPathogenicGenes(exp_data,exp_genes,pathogenic_gene_ranks,pathogenic_genes,sample_i,model,rank_cutoff )
    %find pathogenic genes with rank higher than cut-off
    risky=pathogenic_genes(find(pathogenic_gene_ranks(:,sample_i)>rank_cutoff));
    
    %Set their expression to zero
    [l,risky_genes_index]=ismember(risky, exp_genes);
    risky_genes_index=risky_genes_index(find(risky_genes_index~=0));

    risky_gene_names=exp_genes(risky_genes_index);
    risky_length=length(risky_genes_index);

    before_exp=exp_data(risky_genes_index,sample_i);

    [l,risky_met_genes_index]=ismember(risky_gene_names, model.genes);
    before_exp=before_exp(l);

    risky_met_genes_index=risky_met_genes_index(find(risky_met_genes_index~=0));
    risky_met_genes_length=length(risky_met_genes_index);
    risky_met_genes=model.genes(risky_met_genes_index);
 
    exp_data(risky_genes_index,sample_i)=0;
    newdata=exp_data(:,sample_i);
end



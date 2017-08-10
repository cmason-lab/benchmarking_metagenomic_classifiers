%This is to output what level that we are working on
taxa_operating='genus';

%Directory that sample files are in 
thedir='/Users/gailr/tool-off-test/prill_results/dec_16_16/genus/';
cd(thedir)

%(these must be result of "automated build")
sample_to_compute='Downsample_Deep.100.Mreads_';


    
%Read in BLAST, LMAT, Metaphlan,Diamond, Kraken, and Gottcha results for sample    
    fid=fopen([thedir sample_to_compute 'BlastMeganFiltered' '.txt'],'r');
    blast=textscan(fid, '%d %*s %f %*[^\n]');
    fclose(fid);
    fid=fopen([thedir sample_to_compute 'LMAT' '.txt'],'r');
    lmat=textscan(fid, '%d %*s %f %*[^\n]');
    fclose(fid);
    fid=fopen([thedir sample_to_compute 'Metaphlan' '.txt'],'r');
    metaphlan=textscan(fid, '%d %*s %f %*[^\n]');
    fclose(fid);
    fid=fopen([thedir sample_to_compute 'Gottcha' '.txt'],'r');
    gottcha=textscan(fid, '%d %*s %f %*[^\n]');
    fclose(fid);
    
    fid=fopen([thedir sample_to_compute 'DiamondMegan' '.txt'],'r');
    diamond=textscan(fid, '%d %*s %f %*[^\n]');
    fclose(fid);
    fid=fopen([thedir sample_to_compute 'KrakenFiltered' '.txt'],'r');
    krakenfilt=textscan(fid, '%d %*s %f %*[^\n]');
    fclose(fid);
    
    %put all taxa called together and "join" relative abundances into giant
    %matrix
    temp=[blast{1};lmat{1};metaphlan{1};gottcha{1};diamond{1};krakenfilt{1}];
    %take union first
    alltaxa=unique(temp);
    ra_table=zeros(length(alltaxa),6);
    [temp blast_index alltaxa_index]=intersect(blast{1},alltaxa);
    ra_table(alltaxa_index,1)=blast{2};
    [temp diamond_index alltaxa_index]=intersect(diamond{1},alltaxa);
    ra_table(alltaxa_index,2)=diamond{2};
    [temp gottcha_index alltaxa_index]=intersect(gottcha{1},alltaxa);
    ra_table(alltaxa_index,3)=gottcha{2};
    [temp kraken_index alltaxa_index]=intersect(krakenfilt{1},alltaxa);
    ra_table(alltaxa_index,4)=krakenfilt{2};
    [temp lmat_index alltaxa_index]=intersect(lmat{1},alltaxa);
    ra_table(alltaxa_index,5)=lmat{2};
    [temp metaphlan_index alltaxa_index]=intersect(metaphlan{1},alltaxa);
    ra_table(alltaxa_index,6)=metaphlan{2};
    
    % now we can just use vector math to take the vote
    blastensemble_voting=(ra_table(:,1)>0)+(ra_table(:,3)>0)+(ra_table(:,5)>0)+(ra_table(:,6)>0);
    diamondensemble_voting=(ra_table(:,2)>0)+(ra_table(:,3)>0)+(ra_table(:,4)>0)+(ra_table(:,5)>0)+(ra_table(:,6)>0);
    

    % Do the blast/diamond priority methods
    blast_nowunnormalized_priority=ra_table(:,1);
    blast_nowunnormalized_priority(blast_nowunnormalized_priority<=0)=ra_table(blast_nowunnormalized_priority<=0,5);
    blast_nowunnormalized_priority(blast_nowunnormalized_priority<=0)=ra_table(blast_nowunnormalized_priority<=0,6);
    blast_nowunnormalized_priority(blast_nowunnormalized_priority<=0)=ra_table(blast_nowunnormalized_priority<=0,3);
            
    blast_nowunnormalized_priority(blastensemble_voting<2)=0;
    blast_abundances_priority=blast_nowunnormalized_priority/sum(blast_nowunnormalized_priority);
    
    
    diamond_nowunnormalized_priority=ra_table(:,4);
    diamond_nowunnormalized_priority(diamond_nowunnormalized_priority<=0)=ra_table(diamond_nowunnormalized_priority<=0,5);
    diamond_nowunnormalized_priority(diamond_nowunnormalized_priority<=0)=ra_table(diamond_nowunnormalized_priority<=0,6);
    diamond_nowunnormalized_priority(diamond_nowunnormalized_priority<=0)=ra_table(diamond_nowunnormalized_priority<=0,2);
    diamond_nowunnormalized_priority(diamond_nowunnormalized_priority<=0)=ra_table(diamond_nowunnormalized_priority<=0,3);
            
    diamond_nowunnormalized_priority(diamondensemble_voting<3)=0;
    diamond_abundances_priority=diamond_nowunnormalized_priority/sum(diamond_nowunnormalized_priority);
    
    %output the results
    fid_ss=fopen(char(strcat(sample_to_compute,'blastensemble-priority.txt')),'w');
    fid_di=fopen(char(strcat(sample_to_compute,'diamondensemble-priority.txt')),'w');
        for inc_ss=1:length(alltaxa)
            if blast_abundances_priority(inc_ss)>0
                fprintf(fid_ss,'%d \t %f \t %f \t %s\n',alltaxa(inc_ss), blast_abundances_priority(inc_ss), blast_abundances_priority(inc_ss), taxa_operating);
            end
            if diamond_abundances_priority(inc_ss)>0
                fprintf(fid_di,'%d \t %f \t %f \t %s\n',alltaxa(inc_ss), diamond_abundances_priority(inc_ss), diamond_abundances_priority(inc_ss), taxa_operating);
            end
            end
     fclose(fid_ss);
     fclose(fid_di);
    

    

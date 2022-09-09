function [AA_Clean, AA_Removed, Reads_Clean, Reads_Removed, Reads_Ratio] = Load_NGS_Data_Long_Sequence(FileName, Library, WriteToCSV)
%%%%%%%%
%This script reads all the sequences of a FASTQ file and prepares tables of
%the translated amino acids (AA). In the FASTQ file, there are invalid 
%sequences, and these needs to be removed. Invalid sequences are removed
%based on criteria like a proper 'GGGS' ending and the lack of '*' or 'X' 
%in the translated sequence.
%
%Note that the variable 'Seq_Nucleotides' needs to be changed in the script
%to locate the correct position of the sequence of interest depending on 
%the primers used. This is done in the first segment of the script.
%
%This script was inspired by the work of Brinton et al. "PHASTpep: Analysis 
% Software for Discovery of Cell-Selective Peptides via Phage Display and 
% Next-Generation Sequencing" (https://doi.org/10.1371/journal.pone.0155244)
%
%
%Input:
%   FileName: Name of the FASTQ file to be read as a string.
%
%   Library: PhD library type, either 'PhD7', 'PhD12', or 'PhDC7C'.
%
%   WriteToCSV: Name of the resulting tables to saved as CSV files as a
%   string. If you do not wish to save the data, just use "No" as input.
%
%
%Output:
%   AA_Clean: Cell of three columns with the translated and valid AA in the
%   first column, the number of reads of a particular AA in the second, and
%   their relative contribution in percetange in the third.
%
%   AA_Removed: Similar to AA_Clean, but contains all the invalid reads.
%
%   Reads_Clean: The total number of valid AA reads.
%
%   Reads_Removed: The total number of invalid AA reads.
%
%   Reads_Ratio: The ratio of invalid reads, i.e.,
%   Reads_Removed/(Reads_Clean + Reads_Removed).
%
%
%%%%%%%%


%----- Library chosen decides length of sequences to load -----%
if strcmp(Library, 'PhD7')
    Seq_Length = 7 + 4;
    Seq_Nucleotides = [128:161];
elseif strcmp(Library, 'PhD12')
    Seq_Length = 12 + 4;
    Seq_Nucleotides = [128:176];
elseif strcmp(Library, 'PhDC7C')
    Seq_Length = 2 + 7 + 1 + 4;
    Seq_Nucleotides = [128:170];
else
    error('Use either PhD7, PhD12, or PhDC7C.');
end



%----- Read the sequences from the FASTQ file -----%
disp('Reading in FASTQ file');
[~, Seqs, ~] = fastqread(FileName);



%----- Removing flanking region -----%
disp('Removing flanking regions');
Seqs_Nucleotide = cell(length(Seqs), 1);
Seqs_Error_Length = cell(length(Seqs), 1);

for i = 1 : length(Seqs)
    
    if length(Seqs{i}) >= Seq_Nucleotides(end)
        Seqs_Nucleotide{i} = Seqs{i}(Seq_Nucleotides);
    else
        Seqs_Error_Length{i} = Seqs{i};
    end
    
end

Seqs_Nucleotide = Seqs_Nucleotide(~cellfun('isempty', Seqs_Nucleotide));
Seqs_Error_Length = Seqs_Error_Length(~cellfun('isempty', Seqs_Error_Length));



%----- Convert nucleotides to AA and put into a frequency sorted matrix -----%
disp('Conversion to amino acids');
AA_Array = nt2aa(Seqs_Nucleotide,...
    'AlternativeStartCodons',false,...
    'ACGTonly', false);

disp('Tabulating amino acids');
AA_Cell = tabulate(AA_Array);



%----- Clean up data by ensuring proper end and lack of '*' and 'X' -----%
disp('Cleaning up the sequences.');

%Remove invalid '*' sequences
disp('Removing sequences with *.');
InvIdx1 = strfind(AA_Cell(:, 1), '*');

AA_Removed = AA_Cell(~cellfun('isempty',InvIdx1),:);
AA_Clean = AA_Cell(cellfun('isempty',InvIdx1),:);

%Ensure proper ending, i.e., 'GGGS'
disp('Ensuring proper ending of GGGS.');
InvIdx2 = strfind(AA_Clean(:, 1), 'GGGS');
InvIdx2 = cellfun(@(InvIdx2) find(InvIdx2 == Seq_Length - 3),...
    InvIdx2, 'UniformOutput', false);

AA_Removed = [AA_Removed; AA_Clean(cellfun('isempty', InvIdx2), :)];
AA_Clean = AA_Clean(~cellfun('isempty', InvIdx2), :);

%Remove invalid 'X' sequences
disp('Removing sequences with X.');
InvIdx3 = strfind(AA_Clean(:, 1), 'X');

AA_Removed = [AA_Removed; AA_Clean(~cellfun('isempty', InvIdx3), :)];
AA_Clean = AA_Clean(cellfun('isempty', InvIdx3), :);

%For the 'PhDC7C'library, additional conditions need to be met
if strcmp(Library, 'PhDC7C')
    
    %Need to check for 'AC' beginning
    InvIdx4 = strfind(AA_Clean(:, 1), 'AC');
    InvIdx4 = cellfun(@(InvIdx4) find(InvIdx4 == 1),...
        InvIdx4, 'UniformOutput', false);
    
    AA_Removed = [AA_Removed; AA_Clean(cellfun('isempty', InvIdx4), :)];
    AA_Clean = AA_Clean(~cellfun('isempty', InvIdx4), :);
    
    %and 'C' after sequence
    InvIdx5 = strfind(AA_Clean(:, 1), 'C');
    InvIdx5 = cellfun(@(InvIdx5) find(InvIdx5 == 10),...
        InvIdx5, 'UniformOutput', false);
    
    AA_Removed = [AA_Removed; AA_Clean(cellfun('isempty', InvIdx5), :)];
    AA_Clean = AA_Clean(~cellfun('isempty', InvIdx5), :);
    
end



%----- Sort the cells -----%
AA_Clean = sortrows(AA_Clean, -2);
AA_Removed = sortrows(AA_Removed, -2);

%----- Update relative percentage for clean and removed cells  -----%
AA_Clean(:, 3) = num2cell(100 * cell2mat(AA_Clean(:, 2)) /...
    sum(cell2mat(AA_Clean(:, 2))));
AA_Removed(:, 3) = num2cell(100 * cell2mat(AA_Removed(:, 2)) /...
    sum(cell2mat(AA_Removed(:, 2))));

%----- Calculate the total number of reads -----%
Reads_Clean = sum(cell2mat(AA_Clean(:,2)));
Reads_Removed = sum(cell2mat(AA_Removed(:,2)));
Reads_Ratio = Reads_Removed / (Reads_Clean + Reads_Removed);



%----- Decide upon writing to a CSV file -----%
if ~strcmp(WriteToCSV, 'No')
    
    AA_Clean_Table = table(AA_Clean(:,1),...
        cell2mat(AA_Clean(:,2)),...
        cell2mat(AA_Clean(:,3)),...
        'VariableNames', {'AA', 'Number', 'Relative (%)'});
    AA_Removed_Table = table(AA_Removed(:,1),...
        cell2mat(AA_Removed(:,2)),...
        cell2mat(AA_Removed(:,3)),...
        'VariableNames', {'AA', 'Number', 'Relative (%)'});
    
    writetable(AA_Clean_Table, WriteToCSV + "_AA_Clean.csv");
    writetable(AA_Removed_Table, WriteToCSV + "_AA_Removed.csv");
    
end



end

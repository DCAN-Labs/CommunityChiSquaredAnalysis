function PermuteCrossCorrelationMatrix(m,behavior,test_type,npermutations)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if size(behavior,1) > 1
    number_of_subjects = size(behavior,1);
elseif size(behavior,2) > 1
    number_of_subjects = size(behavior,2);
end

for i = 1:npermutations
    fprintf('%s %i\n','permuting assignments for the groups: permutation #',i);
    permutedir = sprintf('%s','permute',num2str(i),'/');
    totalpermutedir = strcat(pwd,'/',permutedir);
    if isempty(dir(permutedir)) == 0 && size(dir(permutedir),1) > 2
    else                
        command = sprintf('%s','mkdir ',permutedir);
        system(command);
        s = RandStream.create('mt19937ar','seed',sum(100*clock)); %regenerate the random number sequence
        RandStream.setGlobalStream(s);
        R = randperm(number_of_subjects); 
        permuted_behavior = zeros(number_of_subjects,1);
        permuted_behavior(:,1) = behavior(R);
        permuted_corrmat = GenerateCrossCorrelationMatrix(m,permuted_behavior,0,'fisher',test_type);
        permuted_corrmat(:,:,3) = BinarizeMatrix(permuted_corrmat(:,:,2),0.05,'one');
        savematname = strcat('permuted_corrmat',num2str(i),'.mat');
        savebehaviorname = strcat('permuted_behavior',num2str(i),'.mat');
        save(savematname,'permuted_corrmat');
        command = strcat(['mv',' ',savematname,' ',totalpermutedir]);
        system(command);
        save(savebehaviorname,'permuted_behavior');
        command = strcat(['mv',' ',savebehaviorname,' ',totalpermutedir]);
        system(command);
    end
    fclose('all');
end

end


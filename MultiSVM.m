function [result] = MultiSVM(TrainingSet, Label, TestSet)
    u=unique(Label);
    numClasses=length(u);
    result = zeros(length(TestSet(:,1)),1);
    %build models
    for k=1:numClasses
        %Vectorized statement that binarizes Group
        %where 1 is the current class and 0 is all other classes
        TwoGroups = (Label==u(k));
        models(k) = svmtrain(TrainingSet,TwoGroups);
    end
    %classify test cases
    for j=1:size(TestSet,1)
        for k=1:numClasses
            if(svmclassify(models(k),TestSet(j,:)))
                break;
            end
        end
        result(j) = k;
    end
end
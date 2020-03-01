% function [truePositives, falsePositives, falseNegatives] = performance(realMarks, estimatedMarks, eps)
% Calculates TP, FP and FN for estimated marks.
%
% Inputs
% realMarks: Truth source of marks
% estimatedMarks: Marks estimated by an algorithm
% eps: Value at both sides of a mark
%
% Ouputs
% truePositives: Integer describing the number of TP
% falsePositives: Integer describing the number of FP
% falseNegatives: Integer describing the number of FN
function [truePositives, falsePositives, falseNegatives] = performance(realMarks, estimatedMarks, eps)
	truePositives = 0;
	falsePositives = 0;
	falseNegatives = 0;
	
	for i = 1:length(realMarks)
		min = realMarks(i) - eps;
		max = realMarks(i) + eps;
		res = sum((estimatedMarks > min).*(estimatedMarks < max));
		if (res == 0)
			falseNegatives = falseNegatives + 1;
		else
			truePositives = truePositives + 1;
			if (res ~=1)
				falsePositives = falsePositives + 1;
			end
		end
	end

	for i=1:length(estimatedMarks)
		min = estimatedMarks(i) - eps;
		max = estimatedMarks(i) + eps;
		res = sum((realMarks > min).*(realMarks < max));
		if (res == 0)
			falsePositives = falsePositives + 1;
		end
	end
end
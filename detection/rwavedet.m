% function marks = rwavedet(x, threshold, q)
% R wave detector
%
% Inputs
% x: ECG Signal
% threshold: Threshold to detect R waves
% q: Quantity
%
% Outputs
% marks: Indices for R wave marks detected
function marks = rwavedet(x, threshold, q)
  L = length(x);
  marks = zeros(1, L);
  i = 1;
  while i < L
      j = 0;
      max = i;
      while x(i+j) > threshold
          if L > i + j
              if x(i+j) > x(max)
                  max = i + j;
              end
          end
          j = j+1;
      end
      
      if j >= q - 1
        marks(max) = x(max);
      end
      i = i + j + 1;
  end
end

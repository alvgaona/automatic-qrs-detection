% function y = resample(x, N, M, b, a)
% Resampling a signal by N/M ratio.
%
% Inputs
% x: Signal
% N: Upsampling factor
% M: Downsampling factor
% b: FIR filter numerator
% a: FIR filter denominator
%
% Outputs
% y: Resampled signal
function y = resample(x, N, M, b, a)
  % Upsampling
  x_up = zeros(1,length(x)*N);
  for i = 1:length(x)
    x_up(i*N) = x(i);
  end

  x_up = filter(b, a, x_up);

  % Downsampling
  y = zeros(1,floor(length(x_up)/M));
  for i = 1:length(y)
    y(i) = x_up(i*M);
  end
end

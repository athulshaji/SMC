function [output,fsout]=myresample(input,fsin,fsout);

K=gcd(fsin,fsout);
P=fsout/K;
Q=fsin/K;
output=resample(input,P,Q);


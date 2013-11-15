function [ pdadx_out ] = pdadx( xq,pq )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pdadx_out = pq .* dadx( xq );

end


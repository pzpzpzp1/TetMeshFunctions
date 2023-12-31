function qout = paul_quatmultiply( q, varargin )
%  QUATMULTIPLY Calculate the product of two quaternions.
%   N = QUATMULTIPLY( Q, R ) calculates the quaternion product, N, for two
%   given quaternions, Q and R.  Inputs Q and R can be either M-by-4 matrices 
%   containing M quaternions, or a single 1-by-4 quaternion.  N returns an 
%   M-by-4 matrix of quaternion products.  Each element of Q and R must be a
%   real number.  Additionally, Q and R have their scalar number as the first 
%   column.
%
%   Examples:
%
%   Determine the product of two 1-by-4 quaternions:
%      q = [1 0 1 0];
%      r = [1 0.5 0.5 0.75];
%      mult = quatmultiply(q, r)
%
%   Determine the product of a 1-by-4 quaternion with itself:
%      q = [1 0 1 0];
%      mult = quatmultiply(q)
%
%   Determine the product of 1-by-4 and 2-by-4 quaternions:
%      q = [1 0 1 0];
%      r = [1 0.5 0.5 0.75; 2 1 0.1 0.1];
%      mult = quatmultiply(q, r)
%
%   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD, QUATNORM, 
%   QUATNORMALIZE, QUATROTATE.

%   Copyright 2000-2011 The MathWorks, Inc.

%   Note: Quaternion multiplication is not commutative.

narginchk(1, 2);

if any(~isreal(q(:)))
    error(message('aero:quatnorm:isNotReal1'));
end

if (size(q,2) ~= 4)
    error(message('aero:quatnorm:wrongDimension1'));
end

if nargin == 1
    r = q;
else
    r = varargin{1};
    if any(~isreal(r(:)))
        error(message('aero:quatnorm:isNotReal2'));
    end
    if (size(r,1) ~= 4)
        error(message('aero:quatnorm:wrongDimension2'));
    end
    if (size(r,2) ~= size(q,1) && ~( size(r,2) == 1 || size(q,1) == 1))
    %     error(message('aero:quatnorm:wrongDimension3'));
    end
end

% Calculate vector portion of quaternion product
% vec = s1*v2 + s2*v1 + cross(v1,v2)

qout = zeros(size(q,1),size(r,2),4);
qout(:,:,2) = q(:,1).*r(2,:) + r(1,:).*q(:,2) + q(:,3).*r(4,:)-q(:,4).*r(3,:);
qout(:,:,3)= q(:,1).*r(3,:) + r(1,:).*q(:,3) + q(:,4).*r(2,:)-q(:,2).*r(4,:);
qout(:,:,4)= q(:,1).*r(4,:) + r(1,:).*q(:,4) + q(:,2).*r(3,:)-q(:,3).*r(2,:);

%{
vec = [q(:,1).*r(2,:) q(:,1).*r(3,:) q(:,1).*r(4,:)] + ...
         [r(1,:).*q(:,2) r(1,:).*q(:,3) r(1,:).*q(:,4)]+...
         [ q(:,3).*r(4,:)-q(:,4).*r(3,:) ...
           q(:,4).*r(2,:)-q(:,2).*r(4,:) ...
           q(:,2).*r(3,:)-q(:,3).*r(2,:)];
%}

% Calculate scalar portion of quaternion product
% scalar = s1*s2 - dot(v1,v2)
scalar = q(:,1).*r(1,:) - q(:,2).*r(2,:) - ...
             q(:,3).*r(3,:) - q(:,4).*r(4,:);
qout(:,:,1)=scalar;
       

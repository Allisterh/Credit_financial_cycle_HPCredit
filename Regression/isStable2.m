function [indicator,eigenvalues] = isStable2(A)
%ISSTABLE Determine the stability of a lag operator polynomial
%
% Syntax:
%
%	[indicator,eigenvalues] = isStable(A)
%
% Description:
%
%   Given a lag operator polynomial A(L), compute the eigenvalues of the 
%   characteristic polynomial associated with A(L) and determine stability.
%   The stability condition requires that the magnitudes of all eigenvalues 
%   of the characteristic polynomial are less than 1 to within a small
%   numerical tolerance.
%
% Input Argument:
%
%	A - Lag operator polynomial object, as produced by LagOp.
%
% Output Arguments:
%
%   indicator - Boolean indicator for the stability test. TRUE indicates
%     that A(L) is stable and that the magnitudes of all eigenvalues of its
%     characteristic polynomial are less than one (see notes below), and 
%     FALSE indicates otherwise.
%
%	eigenvalues - Eigenvalues of the characteristic polynomial associated
%	  with A(L). The length of eigenvalues is the product of the degree and
%	  dimension of A(L).
%
% Notes:
%
%   o Zero-degree polynomials are always stable.
%
%   o For polynomials of degree greater than zero, the presence of 
%     NaN-valued coefficients returns a FALSE stability indicator and a
%     vector of NaN's in eigenvalues.
%
%   o When testing for stability, the comparison incorporates a small 
%     numerical tolerance. The indicator is TRUE when the magnitudes of all 
%     eigenvalues are less than 1 - 10*eps, where eps is machine precision.
%
%     Users who wish to incorporate their own tolerance (including 0) may
%     simply ignore the indicator and determine stability as follows:
%
%     [~,eigenvalues] = isStable(A);
%     indicator       = all(abs(eigenvalues) < (1 - tol));
%
%     for some small, nonnegative tolerance tol.
%
% Reference:
%
%   [1] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also LagOp/LagOp.

% Copyright 2010 The MathWorks, Inc.

p = A.Degree;            % A(L) polynomial degree

if p == 0                % Zero-degree polynomials are stable
   indicator   = true;
   eigenvalues = [];
   return
end

%
% Form the multivariate companion matrix (F), under the usual assumption
% that the zero-lag coefficient matrix (A0) is an identity. The companion
% matrix is derived by rewriting a degree p polynomial as an equivalent
% polynomial of degree 1 ([1], page 259). 
%

n = A.Dimension;         % Time series dimension
F = zeros(n*p);          % Preallocate the companion matrix

for i = 1:p
    F(1:n,n*i-n+1:n*i) = -(A.Coefficients{i});
end

F((n+1):end,1:(end-n)) = eye(n*(p-1));

%
% Include the effects of a non-identity A0 coefficient by forming a block
% diagonal matrix of the form [A0 I I ...] on the main diagonal and
% computing the generalized eigenvalues of the companion matrix and block
% diagonal matrix. Solving for the generalized eigenvalues of these two
% matrices avoids left-matrix division (\) operations, and improves
% stability of the calculations.
%
% Near singular A0 matrices produce very large eigenvalues which indicate 
% unstable A(L) polynomials, whereas all-zero A0 matrices, which are otherwise 
% stable, produce eigenvalues of infinity which are replaced with zeros.
%
% The following also tests for the presence of NaN-valued coefficients, and
% handles NaN's as a special case returning FALSE as a stability indicator
% and a vector of NaN's as eigenvalues.
%

A0 = A.Coefficients{0};            % Coefficient matrix at lag 0

if any(any(isnan(F))) || any(any(isnan(A0))) || any(any(isinf(F))) || any(any(isinf(A0)))
   indicator   = false;
   eigenvalues = nan(n*p,1);
else
   I  = {eye(n)};                  % Preallocate an identity
   eigenvalues                     = eig(F,blkdiag(A0,I{ones(1,p-1)}));
   eigenvalues(isinf(eigenvalues)) = 0;
   indicator = all(abs(eigenvalues) < (1 - 10*eps));
end
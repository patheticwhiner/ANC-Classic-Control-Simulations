function result = addPolynomials(poly1, poly2, str)
    % ADDPOLYNOMIALS Adds two polynomials and merges like terms.
    %
    %   result = addPolynomials(poly1, poly2)
    %   Adds two polynomials represented by their coefficients and returns
    %   the resulting polynomial with merged like terms.
    %
    %   Inputs:
    %       poly1 - Coefficients of the first polynomial (row vector).
    %       poly2 - Coefficients of the second polynomial (row vector).
    %
    %   str:
    %       left - 在左侧补零
    %       right - 在右侧补零
    %
    %   Output:
    %       result - Coefficients of the resulting polynomial (row vector).

    if nargin < 3
        str = "right";
    end
    % Ensure inputs are row vectors
    poly1 = poly1(:).';
    poly2 = poly2(:).';

    % Determine the lengths of the input polynomials
    len1 = length(poly1);
    len2 = length(poly2);
    % Determine the maximum length
    max_len = max(len1, len2);
    
    % Add the polynomials element-wise
    poly1_padded = padPolynomial(poly1, max_len, str);
    poly2_padded = padPolynomial(poly2, max_len, str);
    result = poly1_padded + poly2_padded;
    
    % Trim the polynomials
    result = trimPolynomial(result, str);
end
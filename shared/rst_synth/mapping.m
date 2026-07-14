function [output1, output2] = mapping(Ao, Bo, So, Ro, input1, input2, mode)
    % 通用函数，根据 mode 的不同，实现两种映射
    % mode = 1: 从 (Ao, Bo, So, Ro, A, B) 计算 (Delta, Gamma)
    % mode = 2: 从 (Ao, Bo, So, Ro, Delta, Gamma) 计算 (A, B)
    if mode == 1
        % 计算 Delta 和 Gamma
        Delta = addPolynomials(conv(input2, Ao), -conv(Bo, input1));
        Gamma = addPolynomials(conv(input1, So), conv(input2, Ro));
        [Delta_coprime, Gamma_coprime] = makeCoprime(Delta, Gamma);
        output1 = Delta_coprime;
        output2 = Gamma_coprime;
    elseif mode == 2
        % 计算 A 和 B
        B = addPolynomials(conv(input2, Bo), conv(input1, So));
        A = addPolynomials(conv(input2, Ao), -conv(input1, Ro));
        output1 = A;
        output2 = B;
    else
        error('Invalid mode. Mode must be 1 or 2.');
    end
end
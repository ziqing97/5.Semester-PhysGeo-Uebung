function Plm = Plm_recursive(l_max, theta, type)
% Legendre functions, recursive,  the index is equal l+1 or m+1
%
%% Initial values
t = cosd(theta);
Plm = cell(l_max+1);
Plm(:, :) = {zeros(size(t))};
%%
if strcmp(type, 'norm')
    % The normalized Plm
    Plm{0+1,0+1} = ones(size(t));
    Plm{1+1,0+1} = sqrt(3)*t;
    Plm{1+1,1+1} = sqrt(3*(1 - t.^2));
    % Recursive
    for l = 2:l_max
        for m = 0:l
            if (m+1) == (l+1)
                Plm{l+1, m+1} = sqrt((2*l+1)/(2*l)) .* sqrt(1-t.^2) .*...
                    Plm{(l-1)+1, (m-1)+1};
            else
                Plm{l+1, m+1} = sqrt((2*l+1)/((l+m)*(l-m))) * ...
                    (sqrt(2*l-1) * t .* Plm{(l-1)+1, m+1} - ...
                    sqrt(((l-1+m)*(l-1-m))/(2*l-3)) .* Plm{(l-2+1), m+1});
            end
        end
    end
else
    % The normal Plm
    Plm{0+1,0+1} = ones(size(t));
    Plm{1+1,0+1} = t;
    Plm{1+1,1+1} = sqrt(1 - t.^2);
    % Recursive
    for l = 2:l_max
        for m = 0:l
            if (m+1) == (l+1)
                Plm{l+1, m+1} = (2*l-1) * sqrt(1 - t.^2) .* Plm{(l-1)+1, (l-1)+1};
            else
                Plm{l+1, m+1} = 1/(l-m) * ((2*l-1) * t.*Plm{(l-1)+1, m+1} - ...
                    (l-1+m) .* Plm{(l-2)+1, m+1});
            end
        end
    end
end

% Change all the zeros into NaN
for i = 1:l_max+1
    for j = i+1:l_max+1
        Plm{i,j} = NaN;
    end
end

end


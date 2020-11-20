function[Plm] = Normalized_Lengendre(l_num,theta)
% mit dieser Funktion kann man normailzed Lendendre rekusiv rechnen
t_num = cosd(theta);
Plm = cell(l_num + 1);
Plm(:,:) = {NaN};
Plm(1,1) = {ones(1,length(theta))};
Plm(2,1) = {sqrt(3) .* t_num};
Plm(2,2) = {sqrt(3 * (1 - t_num.^2))};
for il = 2 : l_num
    for im = 0:il
        if im == il  
            Plm(il + 1, im + 1) = {sqrt((2 * il + 1) ./ (2 * il)) .* sqrt(1 - t_num.^2) .* Plm{il,im}};
        elseif il == im + 1
            Plm(il + 1, im + 1) = {sqrt((2 * il + 1) ./ ((il + im) .* (il - im))) .* sqrt(2 .* il - 1) .* t_num .* Plm{il,im + 1}};
        else
            Plm(il + 1, im + 1) = {sqrt((2 .* il + 1) ./ ((il + im) .* (il - im))) .* (sqrt(2 .* il - 1) .* t_num .* Plm{il,im + 1} - ...
            sqrt((il - 1 + im) .* (il - 1 - im) ./ (2 .* il - 3)) .* Plm{il - 1, im + 1})};
        end
    end
end
end
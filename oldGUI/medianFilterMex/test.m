
image = rand(100);

filter = 5;

[li,lj,lv] = createMedianFilter(size(image),filter);

operator = sparse(li,lj,lv);

result = operator * image(:);
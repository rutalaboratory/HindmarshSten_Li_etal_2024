function prob = dist(array,binsize,low_bound,up_bound);
% for histogram plot
id_dist = low_bound:binsize:up_bound;
prob = [];

for i = 1:length(id_dist)-1
    id_start = id_dist(i);
    id_end = id_dist(i+1);
    
    ids = find(array > id_start & array <= id_end);
    prob(1,i) = length(ids)/length(array);
end

end


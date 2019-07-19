function o = one(w)
w(w<0) = 0;
o = w/sum(w);
end
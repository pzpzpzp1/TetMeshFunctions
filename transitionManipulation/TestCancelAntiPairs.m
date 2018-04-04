transitions = [-1 2 -2 2 3 -3 3 3 3 4 4 -4 1];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == [2 3 3 3 4]));

transitions = [-1];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == [-1]));

transitions = [-1 1];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == []));

transitions = [];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == []));

transitions = [-1 2 1];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == [2]));

transitions = [-1 2 67 1];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == [2 67]));

transitions = [-1 1 2 1];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == [2 1]));

transitions = [-1 1 2 2];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == [2 2]));

transitions = [-1 1];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == []));

transitions = [-5 5 1];
t2 = cancelAntiPairs(transitions);
assert(all(t2 == [1]));

'TestCancelAntiPairs all passed!'
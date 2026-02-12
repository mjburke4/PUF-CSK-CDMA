function code = auth_prf_to_spread(rBits, nonce, L)
% Deterministic PRF-like mapping (simulation only): (rBits, nonce) -> PN code (±1)
% Not cryptographic; fine for simulation/reproducibility.

% Pack bits + nonce to bytes
bytes = [packbits(uint8(rBits)), typecast(uint32(nonce),'uint8')];

% 32-bit FNV-1a hash
h = fnv1a32(bytes);

% Seed a local PRNG and generate ±1 chips
s = RandStream('mt19937ar','Seed', double(h));
code = (randi(s,[0 1], L, 1)*2 - 1);  % ±1 column vector
end

function bytes = packbits(u)
% Pack logical/uint8 vector of bits (0/1) into bytes (LSB first).
u = uint8(u(:)~=0);
n = numel(u); nb = ceil(n/8);
bytes = zeros(1,nb,'uint8');
for i=1:n
    b = bitshift(uint8(1), mod(i-1,8));
    if u(i), bytes(1+floor((i-1)/8)) = bitor(bytes(1+floor((i-1)/8)), b); end
end
end

function h = fnv1a32(bytes)
h = uint32(2166136261); prime = uint32(16777619);
for i=1:numel(bytes)
    h = bitxor(h, uint32(bytes(i)));
    h = uint32( mod( uint64(h) * uint64(prime), 2^32 ) );
end
end


unsigned char getBit(unsigned char w, unsigned char i) {
	return (w >> i) & 0x01;
}

unsigned char setBit(unsigned char w, unsigned char i, unsigned char b) {
	return w | (b<<i);
}
unsigned char xor_all(unsigned char w) {
    unsigned char res = 1;
     for (unsigned char i = 0; i < 8; i++) {
    res = res ^ getBit(w,i);
     }
     return res;
}

unsigned char rule30(unsigned char bl, unsigned char b, unsigned char br) {

unsigned char rule[8]={0,0,0,1,1,1,1,0};

unsigned char idx = (bl<<2)+(b<<1)+br;
return rule[idx];

}


int runCA(int niters) {
    unsigned char w = 188;
    unsigned char wn = 0x0;
    for (int iter=0;iter<niters;iter++) {
        wn = 0x0;
        for (unsigned char i = 1; i < 7; i++) {
            unsigned char im1=i-1;
            unsigned char ip1=i+1;
            unsigned char bl=getBit(w,im1);
            unsigned char br=getBit(w,ip1);
            unsigned char b=getBit(w,i);
            unsigned char bn = rule30(bl,b,br);
            wn = setBit(wn,i,bn);   
        }

        unsigned char bl=getBit(w,7);
        unsigned char br=getBit(w,0);
        unsigned char bn=bl | br;
        wn = setBit(wn,0,bn);  
        wn = setBit(wn,7,bn);  
        w = wn;
            }
    int res = xor_all(w);
    return (res == 0 ? 1 : res);
}

#include <omnetpp.h>
#include <cexception.h>
#include <vector>

#define MAX(a,b) {(a < b)?b:a}
#define MIN(a,b) {(a < b)?a:b}
#define XML_PARSER_ASSERT(expr) if ((expr) == false) {\
        throw cException("Error parsing model description.");\
    }

namespace Utifunc {

    enum process_type {

        process_map
    };

    enum dist_type {

        dist_uniform,

        dist_exponential
    };


    // Matrix class
    class Matrix {
        public:

            // Constructor for empty Matrix
            Matrix();

            // Constructor for Matrix
            // Parameters:
            //  rows    number of rows
            //  cols    number of columns
            Matrix(unsigned int rows, unsigned int cols);

            // Constructor for Matrix
            // Parameters:
            //  m   two-dimensional array with matrix entries
            Matrix(std::vector< std::vector<double> > m);

            // Compute row sum
            double getRowSum(unsigned int row);
            // Get number of rows.
            unsigned int getNumberRows() const;
            // Get number of columns.
            unsigned int getNumberCols() const;
            // Compute inverse matrix.
            Matrix inverse();
            // Transpose matrix.
            Matrix transpose();
            // Access a column of the matrix.
            std::vector<double> getCol(unsigned int col);
            // Replace a column of the matrix.
            void setCol(unsigned int col, std::vector<double> vec);

            // Get a submatrix.
            // Returns the submatrix defined by the rows i1 to i2 and columns j1 to j2
            Matrix getSubmatrix(unsigned int i1, unsigned int i2, unsigned int j1,
                unsigned int j2) const;
            // Insert a submatrix at position i,j.
            void setSubmatrix(unsigned int i, unsigned int j, Matrix B);
            // Access matrix entries.
            std::vector<double>& operator[] (unsigned int i);
            // Access matrix entries.
            double& operator() (unsigned int i, unsigned int j);
            // Matrix Matrix multiplication
            Matrix operator* (Matrix B);
            // Matrix Scalar multiplication
            Matrix operator* (double a);
        private:
            // two-dimensional array with matrix entries
            std::vector< std::vector<double> > m;
    };
    // Matrix subtraction
    Matrix operator- (Matrix A);


        Matrix gauss_jordan(Matrix m);


        Matrix gauss_jordan_inv(Matrix M);
    // ----------------------------------------------------------------------------

    std::pair<Matrix, std::vector<double> > singular_value_decomposition(Matrix M);

    // ----------------------------------------------------------------------------

    // convert double value to string
    std::string toString(double value);



    // abstract super class for distributions
    class Dist {
        public:
            virtual ~Dist();
            // Type of distribution
            dist_type getType();
            // Compute the inverse cumulative distribution function
            virtual double inverse_cdf(double x) = 0;
            // return string representation with description of distribution
            virtual std::string getDescription() = 0;
        protected:
            dist_type type;
    };

    // Uniform distribution
    class Uniform : public Dist {
        public:
            Uniform(double a, double b);
            ~Uniform();
            double inverse_cdf(double x);
            std::string getDescription();
        private:
            // lower bound
            double a;
            // upper bound
            double b;
    };

    // Exponential distribution
    class Exponential : public Dist {
        public:
            Exponential(double rate);
            ~Exponential();
            double inverse_cdf(double x);
            std::string getDescription();
        private:
            // rate parameter
            double rate;
    };



    // Abstract superclass for stochastic processes
    class Process {
        public:
            virtual ~Process();
            // get type of process
            process_type getProcessType();
            // determine hext random variate
            virtual double getNextRandomVariate() = 0;
            // return description of process
            virtual std::string getDescription() = 0;
        protected:
            // type of process
            process_type processtype;
    };

    // Markovian Arrival Process
    class MAP : public Process {
        public:
            MAP(Matrix D0, Matrix D1);
            ~MAP();
            double getNextRandomVariate();
            std::string getDescription();
            // compute stationary vector after an arrival of MAP
            std::vector<double> computeStationaryVector();
        private:
            // matrix D0 with transition rates without an arrival
            Matrix D0;
            // matrix D1 with transition rates with an arrival
            Matrix D1;
            // current state of MAP
            unsigned int state;
            // check if matrices D0 and D1 describe a valid MAP
            bool paramCheck();
    };

};


// Generate arrival events from various stochastic processes
class MAPs : public cSimpleModule {
    public:
        MAPs();
        virtual ~MAPs();
    private:
        cMessage *sendMessageEvent;
        // Transformation of the interarrival times
        cDynamicExpression *transform;
        // Expression for transformation of the interarrival times
        std::string transform_str;
        // Stochastic process
        Utifunc::Process* proc;
        // number of generated arrival events
        unsigned int number;
        // load a MAP from an XML description
        Utifunc::MAP* loadMAP(cXMLElement* xml_map);
    protected:
        virtual void initialize();
        virtual void handleMessage(cMessage *msg);
        virtual void finish();
};



Define_Module(MAPs);

MAPs::MAPs() {
    sendMessageEvent = NULL;
    proc = NULL;
    transform = NULL;
    number = 0;
}

MAPs::~MAPs() {
    cancelAndDelete(sendMessageEvent);
    if (proc != NULL) {delete proc;}
    if (transform != NULL) {delete transform;}
}

void MAPs::initialize() {
    // load arrival process description
    const char* modeltype = par("model").xmlValue()->getTagName();
    if (strcmp(modeltype, "data") != 0) {
        throw cException("Unknown model description for arrival process.");
    }
    if (par("model").xmlValue()->getFirstChildWithTag("map") != NULL) {
        // MAP
        this->proc = this->loadMAP(par("model").xmlValue()->getFirstChildWithTag("map"));
    } else  {
        throw cException("Unknown model description for arrival process.");
    }
    // set display string to show information about the arrival process
    std::stringstream displayItem;
    displayItem << "t=" << this->proc->getDescription();
    getDisplayString().updateWith(displayItem.str().c_str());
    // load transformation function for interarrival times
    if (par("transform").stdstringValue() != "") {
        transform_str = par("transform").stdstringValue();
        transform = new cDynamicExpression();
    }

    // schedule first event
    sendMessageEvent = new cMessage("sendMessageEvent");
    scheduleAt(simTime(), sendMessageEvent);
}

Utifunc::MAP* MAPs::loadMAP(cXMLElement* xml_map) {
    cXMLElement* xml_map_elem = xml_map->getFirstChildWithTag("states");
    XML_PARSER_ASSERT((xml_map_elem != NULL));
    unsigned int map_states = atoi(xml_map_elem->getNodeValue());
    xml_map_elem = xml_map->getFirstChildWithTag("d0");
    XML_PARSER_ASSERT((xml_map_elem != NULL));
    const char* map_d0 = xml_map_elem->getNodeValue();
    xml_map_elem = xml_map->getFirstChildWithTag("d1");
    XML_PARSER_ASSERT((xml_map_elem != NULL));
    const char* map_d1 = xml_map_elem->getNodeValue();
    // construct matrices for MAP
    cStringTokenizer d0_tokenizer(map_d0);
    std::vector<double> d0_elem = d0_tokenizer.asDoubleVector();
    cStringTokenizer d1_tokenizer(map_d1);
    std::vector<double> d1_elem = d1_tokenizer.asDoubleVector();
    if ((d0_elem.size() != map_states * map_states) || (d1_elem.size() != map_states * map_states)) {
        throw cException("Error parsing model description.");
    }
    std::vector<std::vector<double> > D0;
    std::vector<std::vector<double> > D1;
    for (unsigned int i=0; i<map_states; i++) {
        std::vector<double> d0_row(map_states);
        std::vector<double> d1_row(map_states);
        for (unsigned int j=0; j<map_states; j++) {
            d0_row[j] = d0_elem[i*map_states+j];
            d1_row[j] = d1_elem[i*map_states+j];
        }
        D0.push_back(d0_row);
        D1.push_back(d1_row);
    }
    return new Utifunc::MAP(Utifunc::Matrix(D0), Utifunc::Matrix(D1));
}

void MAPs::handleMessage(cMessage *msg) {
    ASSERT(msg==sendMessageEvent);
    // generate arrival
    cMessage *arrival = new cMessage("arrival");
    send(arrival, "out");
    number++;
    // determine next arrival time
    double nextInterarrivalTime = this->proc->getNextRandomVariate();
    // transform interarrival time if transformation function is given
    if (transform != NULL) {
        std::string transform_tmp = transform_str;
        std::size_t transform_pos;
        while ((transform_pos = transform_tmp.find_first_of("$ARRIVAL")) != std::string::npos) {
            // replace with interarrival time
            transform_tmp.replace(transform_pos, 8, Utifunc::toString(nextInterarrivalTime));
        }
        // evaluate expression
        transform->parse(transform_tmp.c_str());
        // get new transformed interarrival time
        nextInterarrivalTime = transform->doubleValue(this);
    }
    // schedule next arrival
    scheduleAt(simTime()+nextInterarrivalTime, sendMessageEvent);
}

void MAPs::finish() {
    EV << "Total jobs generated: " << number << endl;
}



Utifunc::Process::~Process() {
}

Utifunc::process_type Utifunc::Process::getProcessType() {
    return this->processtype;
}

// ----------------------------------------------------------------------------

Utifunc::MAP::MAP(Utifunc::Matrix D0, Utifunc::Matrix D1) {
    this->processtype = Utifunc::process_map;
    this->D0 = D0;
    this->D1 = D1;
    if (!this->paramCheck()) {
        throw cException("Not a valid MAP!");
    }
    // determine initial state according to the stationary distribution of the MAP
    state = 0;
    std::vector<double> pi = this->computeStationaryVector();
    double initStateProb = uniform(0.0,1.0);
    double initStateBound = 0;
    for (unsigned int i=0; i<pi.size(); i++) {
        initStateBound = initStateBound + pi[i];
        if (initStateProb <= initStateBound) {
            state = i;
            break;
        }
    }
}

Utifunc::MAP::~MAP() {
}

bool Utifunc::MAP::paramCheck() {
    // check Matrix-Dimensions
    if (D0.getNumberRows() != D1.getNumberRows()) {return false;}
    if (D0.getNumberRows() != D0.getNumberCols()) {return false;}
    if (D1.getNumberRows() != D1.getNumberCols()) {return false;}
    // check matrix entries
    for (unsigned int i = 0; i<D0.getNumberRows(); i++) {
        double d0ii = 0;
        double sum = 0;
        for (unsigned int j=0; j<D0.getNumberCols(); j++) {
            if (i==j) {
                d0ii = D0[i][j];
                // diagonal elements in D0 have to be negative
                if (d0ii > 0) {return false;}
            } else {
                // off-diagonal elements in D0 have to be non-negative
                if (D0[i][j] < 0) {return false;}
                sum = sum + D0[i][j];
            }
            // elements in D1 have to be non-negative
            if (D1[i][j] < 0) {return false;}
            sum = sum + D1[i][j];
        }
        // diagonal elements have to contain the negative row sums
        // (respecting some rounding error)
        if (fabs(d0ii + sum) > 1e-05) {return false;}
    }
    return true;
}

std::vector<double> Utifunc::MAP::computeStationaryVector() {
    std::vector<double> pi(this->D0.getNumberCols());
    Matrix P = ((-D0).inverse())*D1;
    for (unsigned int i=0; i<P.getNumberRows(); i++) {
        P[i][i] = P[i][i] - 1.0;
    }
    P = P.transpose();
    std::pair<Utifunc::Matrix, std::vector<double> > result = Utifunc::singular_value_decomposition(P);
    double wmax=0.0;
    for(unsigned int j=0; j<result.second.size(); j++) {
        if (result.second[j] > wmax) {
            wmax=result.second[j];
        }
    }
    double epsilon = wmax*1.0e-6;
    // We are looking for entries equal 0 in the vector q (which is returned as result.second).
    // If q[i] == 0, the i-th column in matrix A (which is returned as result.first) is
    // our steady-state vector
    bool found = false;
    for (unsigned int i = 0; i<result.second.size(); i++) {
        if (fabs((result.second[i])) <= epsilon) {
            pi = result.first.getCol(i);
            found = true;
        }
    }
    if (found == false) {
        throw cException("SVD failed.");
    }
    // The steady-state vector has to be normalized such that pi 1 = 1.
    double sum = 0;
    for (unsigned int i = 0; i<pi.size(); i++) {
        sum = pi[i] + sum;
    }
    if (sum != 1) {
        for (unsigned int i = 0; i<pi.size(); i++) {
            pi[i] = pi[i] / sum;
        }
    }
    return pi;
}

double Utifunc::MAP::getNextRandomVariate() {
    double interArrivalTime = 0;
    while (true) {
        // determine next transition of underlying CTMC
        double transitionTime = exponential(-1.0/(this->D0[state][state]));
        interArrivalTime = interArrivalTime + transitionTime;
        double stateProb = uniform(0.0, 1.0);
        double stateBound = 0;
        bool foundnewstate = false;
        for (unsigned int i = 0; i<this->D0.getNumberRows(); i++) {
            if (i != state) {
                stateBound = stateBound + this->D0[state][i] / (-this->D0[state][state]);
                if (stateProb <= stateBound) {
                    // transition according to matrix D0
                    foundnewstate = true;
                    state = i;
                    break;
                }
            }
        }
        if (!foundnewstate) {
            for (unsigned int i = 0; i<this->D1.getNumberRows(); i++) {
                stateBound = stateBound + this->D1[state][i] / (-this->D0[state][state]);
                if (stateProb <= stateBound) {
                    // transition according to matrix D1 => arrival
                    foundnewstate = true;
                    state = i;
                    // exit with next event
                    return(interArrivalTime);
                }
            }
        }
    }
    // should never get here
    return(interArrivalTime);
}

std::string Utifunc::MAP::getDescription() {
    std::stringstream desc;
    desc << "MAP(" << D0.getNumberRows() << ")";
    return desc.str();
}


// ----------------------------------------------------------------------------


#define PARAM_CHECK_INVERSE_CDF if ((x < 0.0) || (x>1.0)) {throw cException("Invalid value for inverse cdf.");}

Utifunc::Dist::~Dist() {
}

Utifunc::dist_type Utifunc::Dist::getType() {
    return this->type;
}

// ----------------------------------------------------------------------------

Utifunc::Uniform::Uniform(double a, double b) {
    if (a > b) {throw cException("Invalid range for uniform distribution");}
    this->a = a;
    this->b = b;
    this->type = Utifunc::dist_uniform;
}

Utifunc::Uniform::~Uniform() {
}

double Utifunc::Uniform::inverse_cdf(double x) {
    PARAM_CHECK_INVERSE_CDF
    return (x*(this->b-this->a) + this->a);
}

std::string Utifunc::Uniform::getDescription() {
    return "Uniform";
}



Utifunc::Matrix::Matrix() {
}

Utifunc::Matrix::Matrix(unsigned int rows, unsigned int cols) {
    // check rows and cols
    if ((rows < 1) || (cols < 1)) {
        throw cException("Invalid matrix dimension.");
    }
    // create Matrix filled with zeros
    for (unsigned int i = 0; i<rows; i++) {
        std::vector<double> temp(cols);
        m.push_back(temp);
    }
}

Utifunc::Matrix::Matrix(std::vector< std::vector<double> > m) {
    // check matrix dimensions
    for (unsigned int i = 0; i<m.size(); i++) {
        if (m[i].size() != m[0].size()) {
            throw cException("Not a matrix.");
        }
    }
    this->m = m;
}

double Utifunc::Matrix::getRowSum(unsigned int row) {
    if (row >= this->getNumberRows()) {
        throw cException("Invalid row or col index.");
    }
    double rowsum = 0.0;
    for (unsigned int i=0; i<this->getNumberCols(); i++) {
        rowsum = rowsum + (*this)(row, i);
    }
    return rowsum;
}

unsigned int Utifunc::Matrix::getNumberRows() const {
    return m.size();
}

unsigned int Utifunc::Matrix::getNumberCols() const {
    if (m.size() == 0) {
        // empty matrix
        return 0;
    }
    return m[0].size();
}


Utifunc::Matrix Utifunc::Matrix::transpose() {
    Matrix B = Matrix(getNumberCols(), getNumberRows());
    // exchange rows and cols
    for (unsigned int i = 0; i<getNumberRows(); i++) {
        for (unsigned int j = 0; j<getNumberCols(); j++) {
            B[j][i] = m[i][j];
        }
    }
    return B;
}

std::vector<double> Utifunc::Matrix::getCol(unsigned int col) {
    if (col >= this->getNumberCols()) {
        throw cException("Invalid index.");
    }
    std::vector<double> colvec;
    // copy values from column to vector
    for (unsigned int i = 0; i<getNumberRows(); i++) {
        colvec.push_back(m[i][col]);
    }
    return colvec;
}

void Utifunc::Matrix::setCol(unsigned int col, std::vector<double> vec) {
    if ((col > this->getNumberCols()) || (vec.size() > this->getNumberRows())) {
        throw cException("Invalid index.");
    }
    // copy values from vector to matrix
    for (unsigned int i = 0; i<getNumberRows(); i++) {
        m[i][col] = vec[i];
    }
}

Utifunc::Matrix Utifunc::Matrix::getSubmatrix(unsigned int i1, unsigned int i2, unsigned int j1,
        unsigned int j2) const {
    // check dimension of submatrix
    if ((i1 >= getNumberRows()) || (i2 >= getNumberRows()) ||
            (j1 >= getNumberCols()) || (j2 >= getNumberCols()) ||
            (j1 > j2) || (i1 > i2)) {
        throw cException("Invalid row or col index.");
    }
    // create submatrix
    Matrix B = Matrix(i2 - i1 + 1, j2 -j1 + 1);
    // copy values from matrix to submatrix
    for (unsigned int i = i1; i<= i2; i++) {
        for (unsigned int j=j1; j<=j2; j++) {
            B[i-i1][j-j1] = m[i][j];
        }
    }
    return B;
}

void Utifunc::Matrix::setSubmatrix(unsigned int i, unsigned int j, Utifunc::Matrix B) {
    // check matrix dimensions
    if ((getNumberRows() < i + B.getNumberRows()) || (getNumberCols() < j + B.getNumberCols())) {
        throw cException("Invalid row or col index.");
    }
    // copy values from submatrix to matrix
    for (unsigned int i1 = i; i1 < i + B.getNumberRows(); i1++) {
        for (unsigned int j1 = j; j1 < j + B.getNumberCols(); j1++) {
            m[i1][j1] = B[i1-i][j1-j];
        }
    }
}

std::vector<double>& Utifunc::Matrix::operator[] (unsigned int i) {
    if (i >= this->getNumberRows() ) {
        throw cException("Invalid matrix index.");
    }
    return m[i];
}

double& Utifunc::Matrix::operator() (unsigned int i, unsigned int j) {
    if ((i >= this->getNumberRows()) || (j >= this->getNumberCols())) {
        throw cException("Invalid row or col index.");
    }
    return m[i][j];
}

Utifunc::Matrix Utifunc::Matrix::operator* (Utifunc::Matrix B) {
    unsigned int mr = this->getNumberRows();
    unsigned int n1 = this->getNumberCols();
    unsigned int n2 = B.getNumberRows();
    unsigned int p = B.getNumberCols();
    // check matrix dimensions
    if (n1 != n2) {
        throw cException("Matrices cannot be multiplied.");
    }
    // perform matrix multiplication
    Matrix C = Matrix(mr, p);
    for (unsigned int i = 0; i<mr; i++) {
        for (unsigned int j=0; j<p; j++) {
            double sum = 0;
            for (unsigned int r = 0; r<n1; r++) {
                sum = sum + m[i][r] * B[r][j];
            }
            C[i][j] = sum;
        }
    }
    return C;
}

Utifunc::Matrix Utifunc::Matrix::operator* (double a) {
    Matrix B = Matrix(getNumberRows(), getNumberCols());
    // multiply all matrix elements with a
    for (unsigned int i = 0; i< getNumberRows(); i++) {
        for (unsigned int j=0; j<getNumberCols(); j++) {
            B[i][j] = a*m[i][j];
        }
    }
    return B;
}

Utifunc::Matrix Utifunc::operator- (Utifunc::Matrix A) {
    Matrix B = A * (-1);
    return B;
}

Utifunc::Matrix Utifunc::Matrix::inverse() {
    return Utifunc::gauss_jordan_inv(*this);
}

Utifunc::Matrix Utifunc::gauss_jordan(Utifunc::Matrix m) {
    // Puts given matrix (2D array) into the Reduced Row Echelon Form.
    // Ported from a Python version written by Jarno Elonen in April 2005, released into Public Domain
    double eps = 1.0e-17;
    unsigned int h = m.getNumberRows();
    unsigned int w = m.getNumberCols();
    Utifunc::Matrix mm = m;
    unsigned int maxrow;
    for (unsigned int y = 0; y < h; y++) {
        maxrow = y;
        for (unsigned int y2=y+1; y2<h; y2++) { // Find max pivot
            if (fabs(mm[y2][y]) > fabs(mm[maxrow][y])) {
                maxrow = y2;
            }
        }
        std::vector<double> temp_row = mm[y];
        mm[y] = mm[maxrow];
        mm[maxrow] = temp_row;
        if (fabs(mm[y][y]) <= eps) { //Singular?
            throw cException("Singular matrix.");
        }
        for (unsigned int y2=y+1; y2<h; y2++) { //Eliminate column y
            double c = mm[y2][y] / mm[y][y];
            for (unsigned int x=y; x<w; x++) {
                mm[y2][x] -= mm[y][x] * c;
            }
        }
    }
    for (int y=h-1; y>=0; y--) { // Backsubstitute
        double c = mm[y][y];
        for (int y2 = 0; y2<y; y2++) {
            for (int x = w-1; x>y-1; x--) {
                mm[y2][x] -=  mm[y][x] * mm[y2][y] / c;
            }
        }
        mm[y][y] /= c;
        for (unsigned int x = h; x < w; x++) { // Normalize row y
            mm[y][x] /= c;
        }
    }
    return mm;
}

Utifunc::Matrix Utifunc::gauss_jordan_inv(Utifunc::Matrix M) {
    if (M.getNumberRows() != M.getNumberCols()) {
        throw cException("Matrix cannot be inverted.");
    }
    Utifunc::Matrix A = Utifunc::Matrix(M.getNumberRows(), 2*M.getNumberCols());
    A.setSubmatrix(0, 0, M);
    // append identity matrix
    Utifunc::Matrix I = Utifunc::Matrix(M.getNumberRows(), M.getNumberCols());
    for (unsigned int i=0; i<M.getNumberRows(); i++) {
        I[i][i] = 1.0;
    }
    A.setSubmatrix(0, M.getNumberCols(), I);
    // Gauss Jordan Elimination
    Utifunc::Matrix B = Utifunc::gauss_jordan(A);
    // return submatrix with inverse of M
    Utifunc::Matrix C = B.getSubmatrix(0, M.getNumberRows()-1, M.getNumberCols(), 2*M.getNumberCols()-1);
    return C;
}

// ----------------------------------------------------------------------------

std::pair<Utifunc::Matrix, std::vector<double> > Utifunc::singular_value_decomposition(Utifunc::Matrix M) {
    unsigned int n = MAX(M.getNumberRows(), M.getNumberCols());
    Matrix A = Matrix(n, n);
    A.setSubmatrix(0, 0, M);
    unsigned int m = n;
    int p = 0;
    double eps = 1.5e-8;
    double tol = 1e-20;
    std::vector<double> q(n);
    int i, j, k, l, l1, n1, np;
    double c, f, g, h, s, x, y, z;
    std::vector<double> e(n);
    // Householder's reduction to bidiagonal form
    g = x = 0; np = n + p;
    for (i = 1; i<=n; i++) {
        e[i-1] = g; s = 0; l = i + 1;
        for (j = i; j<=m; j++) {s = s + A(j-1, i-1) * A(j-1, i-1);}
        if (s < tol) {g = 0;
        } else {
            f = A(i-1, i-1);
            if (f < 0) {g = sqrt(s);} else {g = -sqrt(s);}
            h = f * g - s; A(i-1, i-1) = (f-g);
            for (j = l; j<=np; j++) {
                s = 0;
                for (k = i; k<=m; k++) {s = s + A(k-1, i-1) * A(k-1, j-1);}
                f = s/h;
                for (k = i; k<=m; k++) {A(k-1, j-1) = A(k-1, j-1) + f * A(k-1, i-1);}
            }
        }
        q[i-1] = g; s = 0;
        if (i<=m) {
            for (j = l; j <= n; j++) {s = s + A(i-1, j-1) * A(i-1, j-1);}
        }
        if (s < tol) {g = 0;
        } else {
            f = A(i-1, i);
            if (f < 0) {g = sqrt(s);} else {g = -sqrt(s);}
            h = f*g - s;
            A(i-1, i) = f-g;
            for (j = l; j<=n; j++) {e[j-1] = A(i-1, j-1)/h;}
            for (j = l; j<=m; j++) {
                s = 0;
                for (k = l; k<=n; k++) {s = s + A(j-1, k-1) * A(i-1, k-1);}
                for (k = l; k<=n; k++) {A(j-1, k-1) = A(j-1, k-1) + s * e[k-1];}
            }
        }
        y = fabs(q[i-1]) + fabs(e[i-1]);
        if (y>x) {x = y;}
    }
    // accumulation of right-hand transformations
    for (i = n; i>=1; i--) {
        if (g != 0) {
            h = A(i-1, i) * g;
            for (j = l; j<=n; j++) {A(j-1, i-1) = A(i-1, j-1)/h;}
            for (j = l; j<=n; j++) {
                s = 0;
                for (k = l; k<=n; k++) {s = s + A(i-1, k-1) * A(k-1, j-1);}
                for (k = l; k<=n; k++) {A(k-1, j-1) = A(k-1, j-1) + s * A(k-1, i-1);}
            }
        }
        for (j = l; j<=n; j++) {
            A(i-1, j-1) = 0;
            A(j-1, i-1) = 0;
        }
        A(i-1, i-1) = 1; g = e[i-1]; l = i;
    }
    eps = eps * x; n1 = n+1;
    for (i = m+1; i<=n; i++) {
        for (j = n1; j<=np; j++) {A(i-1, j-1) = 0;}
    }
    // diagonalization to bidiagonal form
    for (k = n; k>=1; k--) {
        test_f_splitting:
        for (l = k; l>=1; l--) {
            if (fabs(e[l-1]) <= eps) {goto test_f_convergence;}
            if (fabs(q[l-2]) <= eps) {goto cancellation;}
        }
        // cancellation of e[l] if l>1
        cancellation:
        c = 0; s = 1; l1 = l-1;
        for (i = l; i<=k; i++) {
            f = s * e[i-1]; e[i-1] = c * e[i-1];
            if (fabs(f) <= eps) {goto test_f_convergence;}
            g = q[i-1]; q[i-1] = h = sqrt(f*f + g*g); c = g/h; s = -f/h;
            for (j = n1; j<=np; j++) {
                y = A(l1-1, j-1); z = A(i-1, j-1);
                A(l1-1, j-1) = c*y + s*z; A(i-1, j-1) = -s*y+c*z;
            }
        }
        test_f_convergence:
        z = q[k-1];
        if (l==k) {goto convergence;}
        // shift from bottom 2x2 minor
        x = q[l-1]; y = q[k-2]; g = e[k-2]; h = e[k-1];
        f = ((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y); g = sqrt(f*f+1);
        if (f<0) {f = ((x-z)*(x+z)+h*(y/(f-g)-h))/x;} else {f = ((x-z)*(x+z)+h*(y/(f+g)-h))/x;}
        // next QR transformation
        c = s = 1;
        for (i = l+1; i<=k; i++) {
            g = e[i-1]; y = q[i-1]; h = s*g; g = c*g;
            e[i-2] = z = sqrt(f*f+h*h); c = f/z; s = h/z;
            f = x*c+g*s; g = -x*s+g*c; h = y*s; y = y*c;
            for (j = 1; j<=n; j++) {
                x = A(j-1, i-2); z = A(j-1, i-1);
                A(j-1, i-2) = x*c+z*s; A(j-1, i-1) = -x*s+z*c;
            }
            q[i-2] = z = sqrt(f*f+h*h); c = f/z; s = h/z;
            f = c*g+s*y; x = -s*g+c*y;
            for (j = n1; j<=np; j++) {
                y = A(i-2, j-1); z = A(i-1, j-1);
                A(i-2, j-1) = c*y+s*z; A(i-1, j-1) = -s*y+c*z;
            }
        }
        e[l-1] = 0; e[k-1] = f; q[k-1] = x;
        goto test_f_splitting;
        convergence:
        if (z<0) {
            // q[k-1] is made non-negative
            q[k-1] = -z;
            for (j=1; j<=n; j++) {A(j-1, k-1) = -(A(j-1, k-1));}
        }
    }
    return (std::pair<Utifunc::Matrix, std::vector<double> > (A, q));
}

// ----------------------------------------------------------------------------

std::string Utifunc::toString(double value) {
    std::stringstream s;
    s << value;
    return s.str();
}




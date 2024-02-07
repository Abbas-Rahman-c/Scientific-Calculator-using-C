#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include<time.h>
#include <stdlib.h>
#include<string.h>

#define SPEED_OF_LIGHT 299792458             // in meters per second
#define GRAVITATIONAL_CONSTANT 6.67430e-11   // in m^3 kg^-1 s^-2
#define PLANCK_CONSTANT 6.62607015e-34       // in m^2 kg / s
#define ELECTRON_MASS 9.10938356e-31         // in kg
#define PROTON_MASS 1.6726219e-27            // in kg
#define NEUTRON_MASS 1.674927471e-27         // in kg
#define ELECTRON_CHARGE -1.602176634e-19     // in Coulombs
#define PERMITTIVITY_OF_VACUUM 8.854187817e-12 // in F/m (Farad per meter)
#define PERMEABILITY_OF_VACUUM 1.2566370614e-6 // in H/m (Henry per meter)
#define BOLTZMANN_CONSTANT 1.380649e-23      // in J/K (Joules per Kelvin)
#define AVOGADRO_NUMBER 6.02214076e23        // in mol^-1
#define GAS_CONSTANT 8.314462618             // in J/(mol·K)
#define STANDARD_GRAVITY 9.80665             // in m/s^2 (acceleration)
#define STANDARD_ATMOSPHERE 101325           // in Pascals
#define STEFAN_BOLTZMANN_CONSTANT 5.670374419e-8 // in W/(m^2·K^4)
#define FARADAY_CONSTANT 96485.33212         // in C/mol (Coulombs per mole)
#define ELEMENTARY_CHARGE 1.602176634e-19    // in Coulombs
#define RYDBERG_CONSTANT 10973731.568160     // in m^-1
#define BOHR_RADIUS 5.29177210903e-11        // in meters
#define BOHR_MAGNETON 9.2740100783e-24       // in J/T
#define FINE_STRUCTURE_CONSTANT 7.2973525693e-3 // dimensionless
#define ALPHA_PARTICLE_MASS 6.644657230e-27  // in kg
#define DEUTERON_MASS 3.343583719e-27        // in kg
#define MUON_MASS 1.883531594e-28            // in kg
#define PI 3.14159265358979323846            // Pi


// Function to calculate the discriminant for quadratic equations
double discr(double a, double b, double c) {
    return b * b - 4 * a * c;//important line
}

// Function to convert degrees to radians
double toRadians(double degree) {
    return degree * (M_PI / 180.0);
}

// Function to calculate the factorial of a number
long long factorial(int n) {
    if (n == 0) return 1;
    long long fact = 1;
    for (int i = 1; i <= n; i++) {
        fact *= i;
    }
    return fact;
}

// Function to calculate the area and circumference of a circle
void calculateCircle() {
    double radius;
    printf("Enter the radius of the circle: ");
    scanf("%lf", &radius);
    printf("Area: %.2lf\n", M_PI * radius * radius);
    printf("Circumference: %.2lf\n", 2 * M_PI * radius);
}

// Function to calculate the area and perimeter of a rectangle
void calculateRectangle() {
    double length, width;
    printf("Enter the length and width of the rectangle: ");
    scanf("%lf %lf", &length, &width);
    printf("Area: %.2lf\n", length * width);
    printf("Perimeter: %.2lf\n", 2 * (length + width));
}

// Function to calculate the area of a triangle
void calculateTriangle() {
    double base, height;
    printf("Enter the base and height of the triangle: ");
    scanf("%lf %lf", &base, &height);
    printf("Area: %.2lf\n", 0.5 * base * height);
}

// Function to calculate the area and perimeter of a square
void calculateSquare() {
    double side;
    printf("Enter the length of the side of the square: ");
    scanf("%lf", &side);
    printf("Area: %.2lf\n", side * side);
    printf("Perimeter: %.2lf\n", 4 * side);
}
double calculateRadius(double x, double y) {
    return sqrt(x * x + y * y);
}
// Function to calculate the derivative of a polynomial
void calculatederivative(int degree, double coefficient[], double result[]) {
    for (int i = 1; i <= degree; i++) {
        result[i - 1] = i * coefficient[i];
    }
    result[degree - 1] = 0;  // The zero-degree term becomes zero
}

// Function to print a polynomial
void printPolynomial(int degree, double coefficient[]) {
    for (int i = degree; i >= 0; i--) {
        if (i == degree) {
            printf("%.2fx^%d", coefficient[i], i);
        } else if (i == 0) {
            if (coefficient[i] < 0) {
                printf(" - %.2f", -coefficient[i]);
            } else {
                printf(" + %.2f", coefficient[i]);
            }
        } else {
            if (coefficient[i] < 0) {
                printf(" - %.2fx^%d", -coefficient[i], i);
            } else {
                printf(" + %.2fx^%d", coefficient[i], i);
            }
        }
    }
    printf("\n");
}// Function to add two vectors
void vectorAddition(double *vector1, double *vector2, double *result, int size) {
    for (int i = 0; i < size; i++) {
        result[i] = vector1[i] + vector2[i];
    }
}

// Function to read a vector from the user
void readVector(double *vector, int size) {
    for (int i = 0; i < size; i++) {
        printf("Enter element %d: ", i + 1);
        scanf("%lf", &vector[i]);
    }
}

// Function to calculate the dot product of two vectors
double dotProduct(double *vector1, double *vector2, int size) {
    double product = 0;
    for (int i = 0; i < size; i++) {
        product += vector1[i] * vector2[i];
    }
    return product;
}
typedef struct {
    double real;
    double imag;
} ComplexNumber;

// Addition of two complex numbers
ComplexNumber addComplex(ComplexNumber a, ComplexNumber b) {
    ComplexNumber result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

// Subtraction of two complex numbers
ComplexNumber subtractComplex(ComplexNumber a, ComplexNumber b) {
    ComplexNumber result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

// Multiplication of two complex numbers
ComplexNumber multiplyComplex(ComplexNumber a, ComplexNumber b) {
    ComplexNumber result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

// Division of two complex numbers
ComplexNumber divideComplex(ComplexNumber a, ComplexNumber b) {
    ComplexNumber result;
    double denominator = b.real * b.real + b.imag * b.imag;
    result.real = (a.real * b.real + a.imag * b.imag) / denominator;
    result.imag = (a.imag * b.real - a.real * b.imag) / denominator;
    return result;
}

// Function to calculate HCF (GCD) of two numbers
int hcf(int a, int b) {
    if (b == 0) {
        return a;
    } else {
        return hcf(b, a % b);
    }
}


// Function to calculate LCM of two numbers
int lcm(int a, int b) {
    return (a / hcf(a, b)) * b;
}
// Function to add two matrices
void addMatrices(double matrix1[3][3], double matrix2[3][3], double result[3][3], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
}

// Function to read a matrix from the user
void readMatrix(double matrix[3][3], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("Enter element [%d][%d]: ", i + 1, j + 1);
            scanf("%lf", &matrix[i][j]);
        }
    }
}

// Function to print a matrix
void printMatrix(double matrix[3][3], int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%.2lf\t", matrix[i][j]);
        }
        printf("\n");
    }
}
void multiplyMatrices(double firstMatrix[3][3], double secondMatrix[3][3], double mult[3][3], int rowFirst, int columnFirst, int rowSecond, int columnSecond) {
    // Initializing elements of matrix mult to 0.
    for (int i = 0; i < rowFirst; ++i) {
        for (int j = 0; j < columnSecond; ++j) {
            mult[i][j] = 0;
        }
    }

    // Multiplying first and second matrices and storing it in mult
    for (int i = 0; i < rowFirst; ++i) {
        for (int j = 0; j < columnSecond; ++j) {
            for (int k = 0; k < columnFirst; ++k) {
                mult[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
            }
        }
    }
}
double determinant(double matrix[3][3], int n) {
    double det = 0;
    double submatrix[3][3];

    if (n == 2)
        return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
    else {
        for (int x = 0; x < n; x++) {
            int subi = 0;
            for (int i = 1; i < n; i++) {
                int subj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == x)
                        continue;
                    submatrix[subi][subj] = matrix[i][j];
                    subj++;
                }
                subi++;
            }
            det += (pow(-1, x) * matrix[0][x] * determinant(submatrix, n - 1));
        }
    }
    return det;
}
void getCofactor(double matrix[3][3], double temp[3][3], int p, int q, int n) {
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            // Copying into temporary matrix only those element which are not in given row and column
            if (row != p && col != q) {
                temp[i][j++] = matrix[row][col];
                // Row is filled, so increase row index and reset col index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

void adjoint(double matrix[3][3], double adj[3][3]) {
    if (3 == 1) {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of matrix[][]
    int sign = 1;
    double temp[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            // Get cofactor of matrix[i][j]
            getCofactor(matrix, temp, i, j, 3);

            // sign of adj[j][i] positive if sum of row and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, 3 - 1));
        }
    }
}
bool inverse(double matrix[3][3], double inverse[3][3]) {
    // Find determinant of matrix[][]
    double det = determinant(matrix, 3);
    if (det == 0) {
        return false;
    }

    // Find adjoint
    double adj[3][3];
    adjoint(matrix, adj);

    // Find Inverse using formula "inverse(matrix) = adjoint(matrix) / determinant(matrix)"
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            inverse[i][j] = adj[i][j] / det;

    return true;
}
// Function to calculate the integral of a polynomial
void calculateIntegral(int degree, double coefficient[], double integralCoefficients[]) {
    for (int i = 0; i <= degree; i++) {
        integralCoefficients[i + 1] = coefficient[i] / (i + 1);
    }
    integralCoefficients[0] = 0; // Constant of integration, can be set to any value
}

// Function to print the integral of a polynomial
void printIntegral(int degree, double integralCoefficients[]) {
    printf("Integral of the polynomial:\n");
    for (int i = degree + 1; i > 0; i--) {
        if (i == degree + 1) {
            printf("%.2fx^%d", integralCoefficients[i], i);
        } else if (i == 1) {
            printf(" + %.2fx", integralCoefficients[i]);
        } else {
            printf(" + %.2fx^%d", integralCoefficients[i], i);
        }
    }
    printf(" + C\n"); // C is the constant of integration
}
#define MAX_QUESTIONS 100
#define MAX_QUESTION_LENGTH 256

typedef struct {
    char question[MAX_QUESTION_LENGTH];
    char answer[MAX_QUESTION_LENGTH];
} QuizQuestion;

QuizQuestion questions[MAX_QUESTIONS];
int totalQuestions = 0;

void loadQuestions(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    while (fscanf(file, "%[^?]?%[^\n]\n", questions[totalQuestions].question, questions[totalQuestions].answer) != EOF) {
        totalQuestions++;
    }

    fclose(file);
}

int askQuestion(QuizQuestion *question) {
    char userAnswer[MAX_QUESTION_LENGTH];
    printf("%s? ", question->question);
    int c;
    while ((c = getchar()) != '\n' && c != EOF) { }
    if (fgets(userAnswer, MAX_QUESTION_LENGTH, stdin) == NULL) {

    }
    userAnswer[strcspn(userAnswer, "\n")] = 0;

    if (strcasecmp(userAnswer, question->answer) == 0) {
        printf("Correct!\n");
        return 1;
    } else {
        printf("Wrong! The correct answer is %s.\n", question->answer);
        return 0;
    }
}

void solveLinearSystem() {
    float a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3;

    printf("Enter coefficients for the first linear equation (a1x + b1y + c1z = d1): ");
    scanf("%f %f %f %f", &a1, &b1, &c1, &d1);

    printf("Enter coefficients for the second linear equation (a2x + b2y + c2z = d2): ");
    scanf("%f %f %f %f", &a2, &b2, &c2, &d2);

    printf("Enter coefficients for the third linear equation (a3x + b3y + c3z = d3): ");
    scanf("%f %f %f %f", &a3, &b3, &c3, &d3);

    // Calculate determinants
    float det = a1 * (b2 * c3 - b3 * c2) - b1 * (a2 * c3- a3 * c2) + c1 * (a2 * b3 - a3 * b2);

    if (det == 0) {
        printf("The system of linear equations has no unique solution.\n");
    } else {
        // Calculate solutions
        float x = (d1 * (b2 * c3 - b3 * c2) - b1 * (d2 * c3 - d3 * c2) + c1 * (d2 * b3 - d3 * b2)) / det;
        float y = (a1 * (d2 * c3 - d3 * c2) - d1 * (a2 * c3 - a3 * c2) + c1 * (a2 * d3 - a3 * d2)) / det;
        float z = (a1 * (b2 * d3 - b3 * d2) - b1 * (a2 * d3 - a3 * d2) + d1 * (a2 * b3 - a3 * b2)) / det;

        // Display the solutions
        printf("The solutions to the system of linear equations are:\n");
        printf("x = %.2f\n", x);
        printf("y = %.2f\n", y);
        printf("z = %.2f\n", z);
    }
}

int main() {
    int choice, i, n, operationType, degree;
    double a, b, c, x, y,result, discriminant;
    bool validInput = false;
    int angleMode = 0;
    char exitoption;

	do{
	printf("Scientific Calculator\n");
	printf("-------------------------------------------------------------------------\n");
	printf("\n");
	
	printf("Choose the type of operation:\n");
	printf("-------------------------------------------------------------------------\n");
	printf("                 |                       |               |                \n");
	printf("1: Arithmetic\t | 2: Trigonometric\t | 3: Geometric\t | 4: Derivation \n");
	printf("                 |                       |               |                \n");
	printf("-------------------------------------------------------------------------\n");
	printf("                 |                       |               |                \n");
	printf("5: Vector Ops\t | 6: Complex Numbers\t | 7: Constants  | 8: Matrices  \n");
	printf("                 |                       |               |                \n");
	printf("-------------------------------------------------------------------------\n");
	printf("                 |                       |                    |            \n");
	printf("9: Inv Trig Ops\t | 10: Hyperbolic Ops    | 11: Integration    | 12: Simultaneous Equation                          \n");
	printf("                 |                       |                    |           \n");
	printf("-------------------------------------------------------------------------\n");
	printf("13: Or you know everything?\n Let's test it.\n");
	printf("-------------------------------------------------------------------------\n");
	printf("Enter your choice : ");
	scanf("%d", &operationType);
	printf("-------------------------------------------------------------------------\n");


    if (operationType == 1) {
        printf("Select an arithmetic operation to perform:\n");
        printf("-------------------------------------------------------------------------\n");
        printf("1: Addition\t2: Subtraction\t3: Multiplication\t4: Division\n");
        printf("-------------------------------------------------------------------------\n");
        printf("5: Square Root\t6: Power\t7: 1/x\t                8: nth Root of x \n");
        printf("-------------------------------------------------------------------------\n");
        printf("9: x^n\t       10: Factorial  \t11: Log(base10)        12: Modulus(int)\n");
        printf("-------------------------------------------------------------------------\n");
		printf("13:Percentage  14: Discriminant 15: Quadratic Roots\n");
		printf("-------------------------------------------------------------------------\n");
		printf("16: HCF        17: LCM\n");
        printf("-------------------------------------------------------------------------\n");
		printf("Enter your choice: ");
        scanf("%d", &choice);
        printf("-------------------------------------------------------------------------\n");

        switch (choice) {
            case 1: // Addition
                printf("Enter the number of elements to add: ");
                scanf("%d", &n);
                result = 0;
                for (i = 0; i < n; i++) {
                    printf("Enter a number: ");
                    scanf("%lf", &x);
                    result += x;
                }
                printf("Result: %.2lf\n", result);
                break;
            case 2: // Subtraction
                printf("Enter the number of elements to subtract: ");
                scanf("%d", &n);
                if (n <= 0) {
                    printf("Invalid number of elements.\n");
                    break;
                }
                printf("Enter the first number: ");
                scanf("%lf", &result);
                for (i = 1; i < n; i++) {
                    printf("Enter a number: ");
                    scanf("%lf", &x);
                    result -= x;
                }
                printf("Result: %.2lf\n", result);
                break;
            case 3: // Multiplication
                printf("Enter the number of elements to multiply: ");
                scanf("%d", &n);
                result = 1;
                for (i = 0; i < n; i++) {
                    printf("Enter a number: ");
                    scanf("%lf", &x);
                    result *= x;
                }
                printf("Result: %.2lf\n", result);
                break;
            case 4: // Division
                printf("Enter the dividend: ");
                scanf("%lf", &result);
                printf("Enter the divisor: ");
                scanf("%lf", &x);
                if (x == 0) {
                    printf("Cannot divide by zero.\n");
                    break;
                }
                result /= x;
                printf("Result: %.2lf\n", result);
                break;
            case 5: // Square Root
                printf("Enter a number: ");
                scanf("%lf", &x);
                if (x < 0) {
                    printf("Cannot compute square root of a negative number.\n");
                    break;
                }
                printf("Result: %.2lf\n", sqrt(x));
                break;
            case 6: // Power
                printf("Enter the base: ");
                scanf("%lf", &x);
                printf("Enter the exponent: ");
                scanf("%lf", &b);
                printf("Result: %.2lf\n", pow(x, b));
                break;
            case 7: // 1/x
                printf("Enter a number: ");
                scanf("%lf", &x);
                if (x == 0) {
                    printf("Cannot divide by zero.\n");
                    break;
                }
                printf("Result: %.2lf\n", 1 / x);
                break;
            case 8: // nth Root of x
                printf("Enter a number: ");
                scanf("%lf", &x);
                printf("Enter the root value (n): ");
                scanf("%lf", &b);
                if (x < 0 && ((int)b % 2 == 0)) {
                    printf("Cannot compute even root of a negative number.\n");
                    break;
                }
                printf("Result: %.2lf\n", pow(x, 1 / b));
                break;
            case 9: // x^n
                printf("Enter a number (x): ");
                scanf("%lf", &x);
                printf("Enter the power (n): ");
                scanf("%lf", &b);
                printf("Result: %.2lf\n", pow(x, b));
                break;
            case 10: // Factorial
                printf("Enter a non-negative integer: ");
                scanf("%d", &n);
                if (n < 0) {
                    printf("Cannot compute factorial of a negative number.\n");
                    break;
                }
                printf("Result: %lld\n", factorial(n));
                break;
            case 11: // Logarithm (base 10)
                printf("Enter a positive number: ");
                scanf("%lf", &x);
                if (x <= 0) {
                    printf("Logarithm undefined for non-positive values.\n");
                    break;
                }
                printf("Result: %.2lf\n", log10(x));
                break;
            case 12: // Modulus (integers only)
                printf("Enter two integers: ");
                scanf("%d %d", &n, &i);
                printf("Result: %d\n", n % i);
                break;
            case 13: // Percentage
                printf("Enter the total value: ");
                scanf("%lf", &x);
                printf("Enter the percentage: ");
                scanf("%lf", &b);
                printf("Result: %.2lf\n", (x * b) / 100);
                break;
            case 14: // Discriminant
                printf("Enter coefficients a, b, c: ");
                scanf("%lf %lf %lf", &a, &b, &c);
                discriminant = discr(a, b, c);
                printf("Result: %.2lf\n", discriminant);
                break;
            case 15: // Quadratic Roots
                printf("Enter coefficients a, b, c: ");
                scanf("%lf %lf %lf", &a, &b, &c);
                discriminant = discr(a, b, c);
                if (discriminant < 0) {
                    printf("Complex roots.\n");
                } else {
                    double root1 = (-b + sqrt(discriminant)) / (2 * a);
                    double root2 = (-b - sqrt(discriminant)) / (2 * a);
                    printf("Roots: %.2lf and %.2lf\n", root1, root2);
                }
                break;
			case 16: // HCF
			    printf("Enter two integers to find HCF: ");
			    int num1, num2;
			    scanf("%d %d", &num1, &num2);
			    printf("HCF: %d\n", hcf(num1, num2));
			    break;
			case 17: // LCM
			    printf("Enter two integers to find LCM: ");
			    scanf("%d %d", &num1, &num2);
			    printf("LCM: %d\n", lcm(num1, num2));
			    break;


            default:
                printf("Invalid operation.\n");
                break;
        }
    } else if (operationType == 2) {
        // Trigonometric operations
        printf("Select a trigonometric operation to perform\n");
        printf("-------------------------------------------------------------------------\n");
        printf("1: Sin(x)\t2: Cos(x)\t3: Tan(x)\t4: Cosec(x)\n");
        printf("-------------------------------------------------------------------------\n");
		printf("5: Sec(x)\t6: Cot(x)\n");
        printf("-------------------------------------------------------------------------\n");
        printf("Enter your choice: ");
        scanf("%d", &choice);
        printf("-------------------------------------------------------------------------\n");
        choice += 15;  // Adjusting for trigonometric cases (16-21)

        printf("Enter the angle in degrees (1) or radians (2): ");
        scanf("%d", &angleMode);
        printf("Enter the angle value: ");
        scanf("%lf", &x);
        validInput = true;

        if (angleMode == 1) {
            x = toRadians(x);  // Convert to radians
        } else if (angleMode != 2) {
            printf("Invalid angle mode.\n");
            validInput = false;
        }

        if (validInput) {
            switch (choice) {
                case 16: // Sin(x)
                    result = sin(x);
                    break;
                case 17: // Cos(x)
                    result = cos(x);
                    break;
                case 18: // Tan(x)
                	if(tan(x)==90){
                		printf("Tan is undefined for this angle.\n");
                		break;
					}
                    result = tan(x);
                    break;
                case 19: // Cosec(x)
                    if (sin(x) == 0) {
                        printf("Cosecant undefined for this angle.\n");
                        break;
                    }
                    result = 1 / sin(x);
                    break;
                case 20: // Sec(x)
                    if (cos(x) == 90) {
                        printf("Secant undefined for this angle.\n");
                        break;
                    }
                    result = 1 / cos(x);
                    break;
                case 21: // Cot(x)
                    if (tan(x) == 0) {
                        printf("Cotangent undefined for this angle.\n");
                        break;
                    }
                    result = 1 / tan(x);
                    break;
                default:
                    printf("Invalid operation.\n");
                    break;
            }
            printf("Result: %.2lf\n", result);
        }
    } else if (operationType == 3) {
        // Geometric operations
        printf("Select a geometric operation to perform\n");
        printf("-------------------------------------------------------------------------\n");
        printf("1: Circle (Area and Circumference)\t2: Rectangle (Area and Perimeter)\n");
        printf("-------------------------------------------------------------------------\n");
		printf("3: Triangle (Area) \t\t\t4: Square (Area and Perimeter)\n");
		printf("-------------------------------------------------------------------------\n");
        printf("5: Radius of circle.\n");
        printf("-------------------------------------------------------------------------\n");
		printf("Enter your choice: ");
        scanf("%d", &choice);
        printf("-------------------------------------------------------------------------\n");
        switch (choice) {
            case 1:
                calculateCircle();
                break;
            case 2:
                calculateRectangle();
                break;
            case 3:
                calculateTriangle();
                break;
            case 4:
                calculateSquare();
                break;
            case 5:
            	 printf("Enter the x-coordinate of the point: ");
			    scanf("%lf", &x);
			    printf("Enter the y-coordinate of the point: ");
			    scanf("%lf", &y);
			
			    double radius = calculateRadius(x, y);
			    printf("The radius of the circle that passes through the point (%.2f, %.2f) is: %.2f\n", x, y, radius);
				break;
            default:
                printf("Invalid choice.\n");
                break;
        }
    } else if (operationType == 4) {
        // Polynomial derivation
        printf("Enter the degree of the polynomial: ");
        scanf("%d", &degree);
        printf("-------------------------------------------------------------------------\n");

        if (degree < 0) {
            printf("Invalid degree. Please enter a non-negative integer.\n");
            return 1; // only for positive degree equations
        }

        // Assigning memory for coefficients array
        double coefficient[degree + 1];

        // Coefficient from the user
        printf("Enter the coefficients in descending order of powers:\n");
        for (int i = degree; i >= 0; i--) {
            printf("Coefficient for x^%d: ", i);
            scanf("%lf", &coefficient[i]);
        }

        printf("\nOriginal polynomial:\n");
        printPolynomial(degree, coefficient);

        double derivativeCoefficients[degree];
        calculatederivative(degree, coefficient, derivativeCoefficients);

        printf("\nDerivative of the polynomial:\n");
        printPolynomial(degree - 1, derivativeCoefficients);
    }else if (operationType == 5) { // Assuming 5 is for Vector Operations
    printf("Select a vector operation to perform\n");
    printf("1: Vector Addition\n");
    printf("2: Dot Product\n");

    printf("Enter your choice: ");
    scanf("%d", &choice);

    int vectorSize;
    printf("Enter the size of the vectors: ");
    scanf("%d", &vectorSize);

    double vector1[vectorSize], vector2[vectorSize], result[vectorSize];

    printf("Enter elements for the first vector:\n");
    readVector(vector1, vectorSize);
    printf("Enter elements for the second vector:\n");
    readVector(vector2, vectorSize);

    switch (choice) {
        case 1: // Vector Addition
            vectorAddition(vector1, vector2, result, vectorSize);
            printf("Resultant Vector: ");
            for (int i = 0; i < vectorSize; i++) {
                printf("%.2lf ", result[i]);
            }
            printf("\n");
            break;
        case 2: // Dot Product
            printf("Dot Product: %.2lf\n", dotProduct(vector1, vector2, vectorSize));
            break;
    
    }
}
else if (operationType == 6) {
    printf("Select a complex number operation to perform\n");
    printf("1: Addition\n");
    printf("2: Subtraction\n");
    printf("3: Multiplication\n");
    printf("4: Division\n");

    printf("Enter your choice: ");
    scanf("%d", &choice);

    ComplexNumber num1, num2, result;

    printf("Enter real and imaginary parts of the first complex number: ");
    scanf("%lf %lf", &num1.real, &num1.imag);
    printf("Enter real and imaginary parts of the second complex number: ");
    scanf("%lf %lf", &num2.real, &num2.imag);

    switch (choice) {
        case 1: // Addition
            result = addComplex(num1, num2);
            break;
        case 2: // Subtraction
            result = subtractComplex(num1, num2);
            break;
        case 3: // Multiplication
            result = multiplyComplex(num1, num2);
            break;
        case 4: // Division
            result = divideComplex(num1, num2);
            break;
        default:
            printf("Invalid operation.\n");
            return 1;
    }

    printf("Result: %.2lf + %.2lfi\n", result.real, result.imag);
}
else if (operationType == 7) {
    printf("Select a constant to view its value:\n");
    printf("1: Speed of Light\n");
    printf("2: Gravitational Constant\n");
    printf("3: Planck Constant\n");
    printf("4: Electron Mass\n");
    printf("5: Proton Mass\n");
    printf("6: Neutron Mass\n");
    printf("7: Electron Charge\n");
    printf("8: Permittivity of Vacuum\n");
    printf("9: Permeability of Vacuum\n");
    printf("10: Boltzmann Constant\n");
    printf("11: Avogadro Number\n");
    printf("12: Gas Constant\n");
    printf("13: Standard Gravity\n");
    printf("14: Standard Atmosphere\n");
    printf("15: Stefan-Boltzmann Constant\n");
    printf("16: Faraday Constant\n");
    printf("17: Elementary Charge\n");
    printf("18: Rydberg Constant\n");
    printf("19: Bohr Radius\n");
    printf("20: Bohr Magneton\n");
    printf("21: Fine Structure Constant\n");
    printf("22: Alpha Particle Mass\n");
    printf("23: Deuteron Mass\n");
    printf("24: Muon Mass\n");
    printf("25: Pi\n");

    int constantChoice;
    printf("Enter your choice: ");
    scanf("%d", &constantChoice);
    switch (constantChoice) {
        case 1:
            printf("Speed of Light: %e m/s\n", SPEED_OF_LIGHT);
            break;
        case 2:
            printf("Gravitational Constant: %e m^3 kg^-1 s^-2\n", GRAVITATIONAL_CONSTANT);
            break;
        case 3:
            printf("Planck Constant: %e m^2 kg / s\n", PLANCK_CONSTANT);
            break;
        case 4:
            printf("Electron Mass: %e kg\n", ELECTRON_MASS);
            break;
        case 5:
            printf("Proton Mass: %e kg\n", PROTON_MASS);
            break;
        case 6:
            printf("Neutron Mass: %e kg\n", NEUTRON_MASS);
            break;
        case 7:
            printf("Electron Charge: %e C\n", ELECTRON_CHARGE);
            break;
        case 8:
            printf("Permittivity of Vacuum: %e F/m\n", PERMITTIVITY_OF_VACUUM);
            break;
        case 9:
            printf("Permeability of Vacuum: %e H/m\n", PERMEABILITY_OF_VACUUM);
            break;
        case 10:
            printf("Boltzmann Constant: %e J/K\n", BOLTZMANN_CONSTANT);
            break;
        case 11:
            printf("Avogadro Number: %e mol^-1\n", AVOGADRO_NUMBER);
            break;
        case 12:
            printf("Gas Constant: %e J/(mol·K)\n", GAS_CONSTANT);
            break;
        case 13:
            printf("Standard Gravity: %e m/s^2\n", STANDARD_GRAVITY);
            break;
        case 14:
            printf("Standard Atmosphere: %e Pa\n", STANDARD_ATMOSPHERE);
            break;
        case 15:
            printf("Stefan-Boltzmann Constant: %e W/(m^2·K^4)\n", STEFAN_BOLTZMANN_CONSTANT);
            break;
        case 16:
            printf("Faraday Constant: %e C/mol\n", FARADAY_CONSTANT);
            break;
        case 17:
            printf("Elementary Charge: %e C\n", ELEMENTARY_CHARGE);
            break;
        case 18:
            printf("Rydberg Constant: %e m^-1\n", RYDBERG_CONSTANT);
            break;
        case 19:
            printf("Bohr Radius: %e m\n", BOHR_RADIUS);
            break;
        case 20:
            printf("Bohr Magneton: %e J/T\n", BOHR_MAGNETON);
            break;
        case 21:
            printf("Fine Structure Constant: %f\n", FINE_STRUCTURE_CONSTANT);
            break;
        case 22:
            printf("Alpha Particle Mass: %e kg\n", ALPHA_PARTICLE_MASS);
            break;
        case 23:
            printf("Deuteron Mass: %e kg\n", DEUTERON_MASS);
            break;
        case 24:
            printf("Muon Mass: %e kg\n", MUON_MASS);  
            break;
        case 25:
            printf("Pi: %f\n", PI);
            break;
        default:
            printf("Invalid choice.\n");
            break;
    }
}	 else if (operationType == 8) { 
        printf("Select a matrix operation to perform\n");
        printf("1: Matrix Addition\n2: Matrix Multiplication\n3: Determinant of a Matrix\n4: Adjoint of a Matrix\n5: Inverse of a Matrix\n");
        
        int matrixChoice;
        printf("Enter your choice: ");
        scanf("%d", &matrixChoice);

        double matrix1[3][3], matrix2[3][3], result[3][3];
        int i, j;

        // Common input for all matrix operations
        printf("Enter elements of the first 3x3 matrix:\n");
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                scanf("%lf", &matrix1[i][j]);

        switch (matrixChoice) {
            case 1: // Matrix Addition
                printf("Enter elements of the second 3x3 matrix:\n");
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        scanf("%lf", &matrix2[i][j]);
                
                addMatrices(matrix1, matrix2, result, 3, 3);
                printMatrix(result, 3, 3);
                break;

            case 2: // Matrix Multiplication
                printf("Enter elements of the second 3x3 matrix:\n");
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        scanf("%lf", &matrix2[i][j]);
                
                multiplyMatrices(matrix1, matrix2, result, 3, 3, 3, 3);
                printMatrix(result, 3, 3);
                break;

            case 3: // Determinant of a Matrix
                printf("Determinant: %.2lf\n", determinant(matrix1, 3));
                break;

            case 4: // Adjoint of a Matrix
                adjoint(matrix1, result);
                printMatrix(result, 3, 3);
                break;

            case 5: // Inverse of a Matrix
                if (inverse(matrix1, result)) {
                    printMatrix(result, 3, 3);
                } else {
                    printf("Inverse doesn't exist (Determinant is 0).\n");
                }
                break;

            default:
                printf("Invalid matrix operation choice.\n");
                break;
        }
    }
else if (operationType == 9) { // Inverse Trigonometric Operations
    printf("Select an inverse trigonometric operation to perform\n");
    printf("1: asin(x)\t2: acos(x)\t3: atan(x)\n");
    
    printf("Enter your choice: ");
    scanf("%d", &choice);
    
    printf("Enter the value: ");
    scanf("%lf", &x);
    
    switch (choice) {
        case 1:
            result = asin(x);
            break;
        case 2:
            result = acos(x);
            break;
        case 3:
            result = atan(x);
            break;
        default:
            printf("Invalid operation.\n");
            break;
    }
    printf("Result: %.2lf\n", result);
} else if (operationType == 10) { // Hyperbolic Operations
    printf("Select a hyperbolic operation to perform\n");
    printf("1: sinh(x)\t2: cosh(x)\t3: tanh(x)\n");
    printf("4: asinh(x)\t5: acosh(x)\t6: atanh(x)\n");
    
    printf("Enter your choice: ");
    scanf("%d", &choice);
    
    printf("Enter the value: ");
    scanf("%lf", &x);
    
    switch (choice) {
        case 1:
            result = sinh(x);
            break;
        case 2:
            result = cosh(x);
            break;
        case 3:
            result = tanh(x);
            break;
        case 4:
            result = asinh(x); 
            break;
        case 5:
            result = acosh(x); 
            break;
        case 6:
            result = atanh(x); 
            break;
        default:
            printf("Invalid operation.\n");
            break;
    }
    printf("Result: %.2lf\n", result);
}

 else if (operationType == 11) {
    printf("Enter the degree of the polynomial: ");
    scanf("%d", &degree);

    double coefficient[degree + 1];
    double integralCoefficients[degree + 2];

    printf("Enter the coefficients in descending order of powers:\n");
    for (int i = degree; i >= 0; i--) {
        printf("Coefficient for x^%d: ", i);
        scanf("%lf", &coefficient[i]);
    }

    calculateIntegral(degree, coefficient, integralCoefficients);
    printIntegral(degree, integralCoefficients);
}
else if(operationType == 12){

    int degree;

    // Ask the user for the degree of the system of equations
    printf("Enter the degree of the system of equations (1 for linear, 2 for two variables, 3 for three variables): ");
    scanf("%d", &degree);

    switch (degree) {
        case 1:
            solveLinearSystem();
            break;
        case 2:
            solveLinearSystem(); // For simplicity, using the same function for a system of two equations
            break;
        case 3:
            solveLinearSystem(); // Using the same function for a system of three equations
            break;
        default:
            printf("Invalid degree. Supported degrees are 1, 2, and 3.\n");
            break;
    }
}
 else if (operationType == 13) {
        // Quiz Game
        loadQuestions("C:\\Users\\humay\\Desktop\\End project\\questions.txt"); // Load questions from a file
        srand(time(NULL)); // Seed for random number generation

        int score = 0;
        for (int i = 0; i < 5 && i < totalQuestions; i++) {
            int questionIndex = rand() % totalQuestions;
            score += askQuestion(&questions[questionIndex]);
        }
	if(score>=4){
		printf("Maybe you do know everything.\n");
	}else{
		printf("Hehe, maybe you need this calculator.\n");
	}
        printf("Your score: %d/5\n", score);   
}
	else {
        printf("Invalid operation type.\n");
   }
        
    printf("Do you want to exit the program.(Y/N)\n");
    scanf(" %c",&exitoption);
    }while(exitoption=='n'||exitoption=='N');
    return 0;
}
   

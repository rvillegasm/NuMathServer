#include "nlohmann/json.hpp"

#include "pistache/endpoint.h"
#include "pistache/router.h"

#include "numath.h"

#include "external/fparser.hh"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Pistache;
using json = nlohmann::json;

std::string globalFunction = "";
std::string globalGFunction = "";
std::string globalSecondGFunction = "";

double function(double x) {
    FunctionParser fp;
    fp.Parse(globalFunction, "x");
    double variables[1] = {x};
    return fp.Eval(variables);
}

double helperFunction(double x) {
    FunctionParser fp;
    fp.Parse(globalGFunction, "x");
    double variables[1] = {x};
    return fp.Eval(variables);
}

double helperFunction2(double x) {
    FunctionParser fp;
    fp.Parse(globalSecondGFunction, "x");
    double variables[1] = {x};
    return fp.Eval(variables);
}

void postIncrementalSearch(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    std::string func = json_body["func"];
    double x0 = json_body["x0"];
    double delta = json_body["delta"];
    int nIter = json_body["nIter"];
    // preprocess info
    globalFunction = func;
    std::vector<std::vector<double>> table;
    // Incr Search logic
    numath::Interval interval = numath::singleVariableEquations::incrementalSearch(function, x0, delta, nIter, table);
    // Prepare the response
    json res;
    if (interval.wasSuccessful) {
        res["first"] = interval.first;
        res["last"] = interval.last;
    }
    else {
        res["error"] = "Could not find anything";
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postBisection(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    std::string func = json_body["func"];
    double xi = json_body["xi"];
    double xu = json_body["xu"];
    int nIter = json_body["nIter"];
    double tol = json_body["tol"];
    std::string err = json_body["err"];
    // preprocess info
    globalFunction = func;
    std::vector<std::vector<double>> table;
    err = (err == "Absolute") ? "abs" : "rel";
    // Prepare the response
    json res;
    // bisection logic
    try {
        double root = numath::singleVariableEquations::bisection(function, xi, xu, nIter, tol, err.c_str(), table);
        res["root"] = root;
    }
    catch (numath::IntervalException &intEx) {
        res["error"] = intEx.what();
    }
    catch (numath::IterException &iterEx) {
        res["error"] = iterEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postFalsePosition(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    std::string func = json_body["func"];
    double xi = json_body["xi"];
    double xu = json_body["xu"];
    int nIter = json_body["nIter"];
    double tol = json_body["tol"];
    std::string err = json_body["err"];
    // preprocess info
    globalFunction = func;
    std::vector<std::vector<double>> table;
    err = (err == "Absolute") ? "abs" : "rel";
    // Prepare the response
    json res;
    // bisection logic
    try {
        double root = numath::singleVariableEquations::falsePosition(function, xi, xu, nIter, tol, err.c_str(), table);
        res["root"] = root;
    }
    catch (numath::IntervalException &intEx) {
        res["error"] = intEx.what();
    }
    catch (numath::IterException &iterEx) {
        res["error"] = iterEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postFixedPoint(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    std::string func = json_body["func"];
    std::string gFunc = json_body["gFunc"];
    double xa = json_body["xa"];
    int nIter = json_body["nIter"];
    double tol = json_body["tol"];
    std::string err = json_body["err"];
    // preprocess info
    globalFunction = func;
    globalGFunction = gFunc;

    std::vector<std::vector<double>> table;
    err = (err == "Absolute") ? "abs" : "rel";
    // Prepare the response
    json res;
    // bisection logic
    try {
        double root = numath::singleVariableEquations::fixedPoint(function, helperFunction, xa, nIter, tol, err.c_str(), table);
        res["root"] = root;
    }
    catch (numath::IterException &iterEx) {
        res["error"] = iterEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postNewtonSingle(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    std::string func = json_body["func"];
    std::string dFunc = json_body["dFunc"];
    double x0 = json_body["x0"];
    int nIter = json_body["nIter"];
    double tol = json_body["tol"];
    std::string err = json_body["err"];
    // preprocess info
    globalFunction = func;
    globalGFunction = dFunc;

    std::vector<std::vector<double>> table;
    err = (err == "Absolute") ? "abs" : "rel";
    // Prepare the response
    json res;
    // bisection logic
    try {
        double root = numath::singleVariableEquations::newton(function, helperFunction, x0, nIter, tol, err.c_str(), table);
        res["root"] = root;
    }
    catch (numath::IterException &iterEx) {
        res["error"] = iterEx.what();
    }
    catch (numath::DerivativeException &derivEx) {
        res["error"] = derivEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postSecant(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    std::string func = json_body["func"];
    double x0 = json_body["x0"];
    double x1 = json_body["x1"];
    int nIter = json_body["nIter"];
    double tol = json_body["tol"];
    std::string err = json_body["err"];
    // preprocess info
    globalFunction = func;

    std::vector<std::vector<double>> table;
    err = (err == "Absolute") ? "abs" : "rel";
    // Prepare the response
    json res;
    // bisection logic
    try {
        double root = numath::singleVariableEquations::secant(function, x0, x1, nIter, tol, err.c_str(), table);
        res["root"] = root;
    }
    catch (numath::IterException &iterEx) {
        res["error"] = iterEx.what();
    }
    catch (numath::DerivativeException &derivEx) {
        res["error"] = derivEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postMultipleRoots(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    std::string func = json_body["func"];
    std::string dFunc = json_body["dFunc"];
    std::string d2Func = json_body["d2Func"];
    double x0 = json_body["x0"];
    int nIter = json_body["nIter"];
    double tol = json_body["tol"];
    std::string err = json_body["err"];
    // preprocess info
    globalFunction = func;
    globalGFunction = dFunc;
    globalSecondGFunction = d2Func;

    std::vector<std::vector<double>> table;
    err = (err == "Absolute") ? "abs" : "rel";
    // Prepare the response
    json res;
    // bisection logic
    try {
        double root = numath::singleVariableEquations::multipleRoots(function, helperFunction, helperFunction2, x0, nIter, tol, err.c_str(), table);
        res["root"] = root;
    }
    catch (numath::IterException &iterEx) {
        res["error"] = iterEx.what();
    }
    catch (numath::DenominatorException &denomEx) {
        res["error"] = denomEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postSimpleGaussianElimination(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numEq = json_body["numEq"];
    std::vector<double> nums = json_body["nums"];
    // preprocess info
    std::vector<std::vector<double>> augmentedMatrix;
    for (unsigned int i = 0; i < numEq; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < numEq + 1; j++) {
            row.push_back(nums[i*(numEq+1)+j]);
        }
        augmentedMatrix.push_back(row);
    }
    // Prepare the response
    json res;
    // bisection logic
    try {
        std::vector<double> results = numath::systemsOfEquations::simpleGaussianElimination(augmentedMatrix);
        json j_vec(results);
        res["results"] = results;
    }
    catch (numath::DenominatorException &denomEx) {
        res["error"] = denomEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postGaussianEliminationPartialPivot(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numEq = json_body["numEq"];
    std::vector<double> nums = json_body["nums"];
    // preprocess info
    std::vector<std::vector<double>> augmentedMatrix;
    for (unsigned int i = 0; i < numEq; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < numEq + 1; j++) {
            row.push_back(nums[i*(numEq+1)+j]);
        }
        augmentedMatrix.push_back(row);
    }

    // Prepare the response
    json res;
    // bisection logic
    try {
        std::vector<double> results = numath::systemsOfEquations::gaussianEliminationPartialPivot(augmentedMatrix);
        json j_vec(results);
        res["results"] = results;
    }
    catch (numath::DenominatorException &denomEx) {
        res["error"] = denomEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postGaussianEliminationTotalPivot(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numEq = json_body["numEq"];
    std::vector<double> nums = json_body["nums"];
    // preprocess info
    std::vector<std::vector<double>> augmentedMatrix;
    for (unsigned int i = 0; i < numEq; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < numEq + 1; j++) {
            row.push_back(nums[i*(numEq+1)+j]);
        }
        augmentedMatrix.push_back(row);
    }

    // Prepare the response
    json res;
    // bisection logic
    try {
        std::vector<double> results = numath::systemsOfEquations::gaussianEliminationTotalPivot(augmentedMatrix);
        json j_vec(results);
        res["results"] = results;
    }
    catch (numath::DenominatorException &denomEx) {
        res["error"] = denomEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postDoolittle(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numEq = json_body["numEq"];
    std::vector<double> numsA = json_body["numsA"];
    std::vector<double> numsB = json_body["numsB"];
    // preprocess info
    std::vector<std::vector<double>> matrix;
    for (unsigned int i = 0; i < numEq; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < numEq; j++) {
            row.push_back(numsA[i*(numEq)+j]);
        }
        matrix.push_back(row);
    }

    // Prepare the response
    json res;
    // bisection logic
    try {
        std::vector<double> results = numath::systemsOfEquations::doolittleMethod(matrix, numsB);
        json j_vec(results);
        res["results"] = results;
    }
    catch (numath::DenominatorException &denomEx) {
        res["error"] = denomEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postCrout(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numEq = json_body["numEq"];
    std::vector<double> numsA = json_body["numsA"];
    std::vector<double> numsB = json_body["numsB"];
    // preprocess info
    std::vector<std::vector<double>> matrix;
    for (unsigned int i = 0; i < numEq; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < numEq; j++) {
            row.push_back(numsA[i*(numEq)+j]);
        }
        matrix.push_back(row);
    }

    // Prepare the response
    json res;
    // bisection logic
    try {
        std::vector<double> results = numath::systemsOfEquations::croutMethod(matrix, numsB);
        json j_vec(results);
        res["results"] = results;
    }
    catch (numath::DenominatorException &denomEx) {
        res["error"] = denomEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postCholesky(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numEq = json_body["numEq"];
    std::vector<double> numsA = json_body["numsA"];
    std::vector<double> numsB = json_body["numsB"];
    // preprocess info
    std::vector<std::vector<double>> matrix;
    for (unsigned int i = 0; i < numEq; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < numEq; j++) {
            row.push_back(numsA[i*(numEq)+j]);
        }
        matrix.push_back(row);
    }

    // Prepare the response
    json res;
    // bisection logic
    try {
        std::vector<double> results = numath::systemsOfEquations::choleskyMethod(matrix, numsB);
        json j_vec(results);
        res["results"] = results;
    }
    catch (numath::DenominatorException &denomEx) {
        res["error"] = denomEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postGaussSeidel(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numEq = json_body["numEq"];
    double nIter = json_body["nIter"];
    double tol = json_body["tol"];
    double lambda = json_body["lambda"];
    std::string err = json_body["err"];
    std::vector<double> initVals = json_body["initVals"];
    std::vector<double> numsA = json_body["numsA"];
    std::vector<double> numsB = json_body["numsB"];
    // preprocess info
    std::vector<std::vector<double>> table;
    std::vector<std::vector<double>> matrix;
    for (unsigned int i = 0; i < numEq; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < numEq; j++) {
            row.push_back(numsA[i*(numEq)+j]);
        }
        matrix.push_back(row);
    }

    // Prepare the response
    json res;
    // bisection logic
    try {
        std::vector<double> results = numath::systemsOfEquations::solveIterative(matrix, numsB, initVals, tol, nIter, numath::systemsOfEquations::gaussSeidel, err == "Absolute"?numath::absNorm:numath::relNorm, lambda, table);
        json j_vec(results);
        res["results"] = results;
    }
    catch (numath::IterException &iterEx) {
        res["error"] = iterEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postJacobi(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numEq = json_body["numEq"];
    double nIter = json_body["nIter"];
    double tol = json_body["tol"];
    double lambda = json_body["lambda"];
    std::string err = json_body["err"];
    std::vector<double> initVals = json_body["initVals"];
    std::vector<double> numsA = json_body["numsA"];
    std::vector<double> numsB = json_body["numsB"];
    // preprocess info
    std::vector<std::vector<double>> table;
    std::vector<std::vector<double>> matrix;
    for (unsigned int i = 0; i < numEq; i++) {
        std::vector<double> row;
        for (unsigned int j = 0; j < numEq; j++) {
            row.push_back(numsA[i*(numEq)+j]);
        }
        matrix.push_back(row);
    }

    // Prepare the response
    json res;
    // bisection logic
    try {
        std::vector<double> results = numath::systemsOfEquations::solveIterative(matrix, numsB, initVals, tol, nIter, numath::systemsOfEquations::jacobi, err == "Absolute"?numath::absNorm:numath::relNorm, lambda, table);
        json j_vec(results);
        res["results"] = results;
    }
    catch (numath::IterException &iterEx) {
        res["error"] = iterEx.what();
    }
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

void postLagrange(const Rest::Request &request, Http::ResponseWriter response) {
    std::string body = request.body();
    // parse json
    auto json_body = json::parse(body);
    double numPoints = json_body["numPoints"];
    std::vector<double> points = json_body["points"];
    // preprocess info
    std::vector<numath::Point> pointsVec;
    for (unsigned int i = 0; i < numPoints; i+=2) {
        numath::Point p = {points[i], points[i+1]};
        pointsVec.push_back(p);
    }
    // Prepare the response
    json res;
    // bisection logic
    std::string pol = numath::interpolation::lagrange(pointsVec);
    res["pol"] = pol;
    // Send the response
    auto mime = MIME(Application, Json);
    response.send(Http::Code::Ok, res.dump(), mime);
}

int main() {
    Address addr(Ipv4::any(), Port(9080));

    auto opts = Http::Endpoint::options().threads(1);
    Http::Endpoint server(addr);
    server.init(opts);

    Rest::Router router;
    Rest::Routes::Post(router, "/methods/incrSearch", Rest::Routes::bind(postIncrementalSearch));
    Rest::Routes::Post(router, "/methods/bisection", Rest::Routes::bind(postBisection));
    Rest::Routes::Post(router, "/methods/falsePosition", Rest::Routes::bind(postFalsePosition));
    Rest::Routes::Post(router, "/methods/fixedPoint", Rest::Routes::bind(postFixedPoint));
    Rest::Routes::Post(router, "/methods/newtonSingle", Rest::Routes::bind(postNewtonSingle));
    Rest::Routes::Post(router, "/methods/secant", Rest::Routes::bind(postSecant));
    Rest::Routes::Post(router, "/methods/multipleRoots", Rest::Routes::bind(postMultipleRoots));
    Rest::Routes::Post(router, "/methods/SGE", Rest::Routes::bind(postSimpleGaussianElimination));
    Rest::Routes::Post(router, "/methods/GEP", Rest::Routes::bind(postGaussianEliminationPartialPivot));
    Rest::Routes::Post(router, "/methods/GET", Rest::Routes::bind(postGaussianEliminationTotalPivot));
    Rest::Routes::Post(router, "/methods/doolittle", Rest::Routes::bind(postDoolittle));
    Rest::Routes::Post(router, "/methods/crout", Rest::Routes::bind(postCrout));
    Rest::Routes::Post(router, "/methods/cholesky", Rest::Routes::bind(postCholesky));
    Rest::Routes::Post(router, "/methods/gaussSeidel", Rest::Routes::bind(postGaussSeidel));
    Rest::Routes::Post(router, "/methods/jacobi", Rest::Routes::bind(postJacobi));
    Rest::Routes::Post(router, "/methods/lagrange", Rest::Routes::bind(postLagrange));

    server.setHandler(router.handler());
    // server.setHandler(std::make_shared<helloHandler>());
    server.serve();

    server.shutdown();
    // Http::listenAndServe<HelloHandler>("*:9080");
}
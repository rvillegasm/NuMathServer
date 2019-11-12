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

double function(double x) {
    FunctionParser fp;
    fp.Parse(globalFunction, "x");
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
    res["first"] = interval.first;
    res["last"] = interval.last;

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
    Rest::Routes::Post(router, "/methods", Rest::Routes::bind(postIncrementalSearch));

    server.setHandler(router.handler());
    // server.setHandler(std::make_shared<helloHandler>());
    server.serve();

    server.shutdown();
    // Http::listenAndServe<HelloHandler>("*:9080");
}
#include "preprocessor.hpp"

int main(int argc, char* argv[])
{
    const auto scope_guard = lstr::L3sterScopeGuard{argc, argv};
    const auto comm        = std::make_shared< lstr::MpiComm >(MPI_COMM_WORLD);

    Preprocessor preprocessor(comm);
}
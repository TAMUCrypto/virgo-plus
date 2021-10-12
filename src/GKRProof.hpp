//
// Created by 69029 on 8/23/2021.
//

#ifndef VIRGO_PLUS_GKRPROOF_HPP
#define VIRGO_PLUS_GKRPROOF_HPP

using std::istream;

struct GKRProof {
    vector<F> final_claims_u;
    vector<vector<F>> final_claims_v;
    vector<F> final_claims;
    vector<vector<quadratic_poly>> polys_u;
    vector<vector<quadratic_poly>> polys_v;
    vector<vector<quadratic_poly>> polys;

#ifdef USE_VIRGO
    // the proof to open the commitment
    virgo::poly_commit::PolyProof poly_proof;
#endif

    void write(ostream &out) const {
        size_t p0 = out.tellp();
        u64 len = final_claims_u.size();
        out.write(reinterpret_cast<const char *>(&len), sizeof(len));
        out.write(reinterpret_cast<const char *>(final_claims_u.data()), final_claims_u.size() * sizeof(final_claims_u[0]));
        len = final_claims.size();
        out.write(reinterpret_cast<const char *>(&len), sizeof(len));
        out.write(reinterpret_cast<const char *>(final_claims.data()), final_claims.size() * sizeof(final_claims[0]));
        len = final_claims_v.size();
        out.write(reinterpret_cast<const char *>(&len), sizeof(len));
        for (const auto &claim_v: final_claims_v) {
            len = claim_v.size();
            out.write(reinterpret_cast<const char *>(&len), sizeof(len));
            out.write(reinterpret_cast<const char *>(claim_v.data()), claim_v.size() * sizeof(claim_v[0]));
        }
        len = polys_u.size();
        out.write(reinterpret_cast<const char *>(&len), sizeof(len));
        for (const auto &poly: polys_u) {
            len = poly.size();
            out.write(reinterpret_cast<const char *>(&len), sizeof(len));
            out.write(reinterpret_cast<const char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        len = polys_v.size();
        out.write(reinterpret_cast<const char *>(&len), sizeof(len));
        for (const auto &poly: polys_v) {
            len = poly.size();
            out.write(reinterpret_cast<const char *>(&len), sizeof(len));
            out.write(reinterpret_cast<const char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        len = polys.size();
        out.write(reinterpret_cast<const char *>(&len), sizeof(len));
        for (const auto &poly: polys) {
            len = poly.size();
            out.write(reinterpret_cast<const char *>(&len), sizeof(len));
            out.write(reinterpret_cast<const char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        size_t p1 = out.tellp();
        fprintf(stderr, "ps (gkr) = %llu (KB)\n", (p1 - p0) >> 10);
        poly_proof.write(out);
        size_t p2 = out.tellp();
        fprintf(stderr, "ps (pc) = %llu (KB)\n", (p2 - p1) >> 10);
    }

    void send(NetIO *io) {
        u64 len = final_claims_u.size();
        io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
        io->send_data(reinterpret_cast<const char *>(final_claims_u.data()), final_claims_u.size() * sizeof(final_claims_u[0]));
        len = final_claims.size();
        io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
        io->send_data(reinterpret_cast<const char *>(final_claims.data()), final_claims.size() * sizeof(final_claims[0]));
        len = final_claims_v.size();
        io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
        for (const auto &claim_v: final_claims_v) {
            len = claim_v.size();
            io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
            io->send_data(reinterpret_cast<const char *>(claim_v.data()), claim_v.size() * sizeof(claim_v[0]));
        }
        len = polys_u.size();
        io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
        for (const auto &poly: polys_u) {
            len = poly.size();
            io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
            io->send_data(reinterpret_cast<const char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        len = polys_v.size();
        io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
        for (const auto &poly: polys_v) {
            len = poly.size();
            io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
            io->send_data(reinterpret_cast<const char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        len = polys.size();
        io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
        for (const auto &poly: polys) {
            len = poly.size();
            io->send_data(reinterpret_cast<const char *>(&len), sizeof(len));
            io->send_data(reinterpret_cast<const char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        poly_proof.send(io);
    }

    void read(istream &in) {
        size_t p0 = in.tellg();
        u64 len;
        in.read(reinterpret_cast<char *>(&len), sizeof(len));
        final_claims_u.resize(len);
        in.read(reinterpret_cast<char *>(final_claims_u.data()), final_claims_u.size() * sizeof(final_claims_u[0]));
        in.read(reinterpret_cast<char *>(&len), sizeof(len));
        final_claims.resize(len);
        in.read(reinterpret_cast<char *>(final_claims.data()), final_claims.size() * sizeof(final_claims[0]));
        in.read(reinterpret_cast<char *>(&len), sizeof(len));
        final_claims_v.resize(len);
        for (auto &claim_v: final_claims_v) {
            in.read(reinterpret_cast<char *>(&len), sizeof(len));
            claim_v.resize(len);
            in.read(reinterpret_cast<char *>(claim_v.data()), claim_v.size() * sizeof(claim_v[0]));
        }
        in.read(reinterpret_cast<char *>(&len), sizeof(len));
        polys_u.resize(len);
        for (auto &poly: polys_u) {
            in.read(reinterpret_cast<char *>(&len), sizeof(len));
            poly.resize(len);
            in.read(reinterpret_cast<char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        in.read(reinterpret_cast<char *>(&len), sizeof(len));
        polys_v.resize(len);
        for (auto &poly: polys_v) {
            in.read(reinterpret_cast<char *>(&len), sizeof(len));
            poly.resize(len);
            in.read(reinterpret_cast<char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        in.read(reinterpret_cast<char *>(&len), sizeof(len));
        polys.resize(len);
        for (auto &poly: polys) {
            in.read(reinterpret_cast<char *>(&len), sizeof(len));
            poly.resize(len);
            in.read(reinterpret_cast<char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        size_t p1 = in.tellg();
        fprintf(stderr, "ps (gkr) = %llu (KB)\n", (p1 - p0) >> 10);
        poly_proof.read(in);
        size_t p2 = in.tellg();
        fprintf(stderr, "ps (pc) = %llu (KB)\n", (p2 - p1) >> 10);
    }

    void recv(NetIO *io) {
        u64 len;
        io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
        final_claims_u.resize(len);
        io->recv_data(reinterpret_cast<char *>(final_claims_u.data()), final_claims_u.size() * sizeof(final_claims_u[0]));
        io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
        final_claims.resize(len);
        io->recv_data(reinterpret_cast<char *>(final_claims.data()), final_claims.size() * sizeof(final_claims[0]));
        io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
        final_claims_v.resize(len);
        for (auto &claim_v: final_claims_v) {
            io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
            claim_v.resize(len);
            io->recv_data(reinterpret_cast<char *>(claim_v.data()), claim_v.size() * sizeof(claim_v[0]));
        }
        io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
        polys_u.resize(len);
        for (auto &poly: polys_u) {
            io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
            poly.resize(len);
            io->recv_data(reinterpret_cast<char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
        polys_v.resize(len);
        for (auto &poly: polys_v) {
            io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
            poly.resize(len);
            io->recv_data(reinterpret_cast<char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
        polys.resize(len);
        for (auto &poly: polys) {
            io->recv_data(reinterpret_cast<char *>(&len), sizeof(len));
            poly.resize(len);
            io->recv_data(reinterpret_cast<char *>(poly.data()), poly.size() * sizeof(poly[0]));
        }
        poly_proof.recv(io);
    }
};

#endif //VIRGO_PLUS_GKRPROOF_HPP

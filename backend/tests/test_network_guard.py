"""The test-wide network guard (tests/conftest.py) must actually refuse
outbound connections -- otherwise it reads as protection while CI keeps
downloading every real cohort -- and must not break loopback."""
import socket
import threading

import httpx
import pytest

from tests.conftest import NetworkBlockedInTests


def test_outbound_http_is_refused_fast():
    with pytest.raises((httpx.ConnectError, NetworkBlockedInTests)):
        httpx.get("https://ftp.ncbi.nlm.nih.gov/geo/", timeout=30)


def test_raw_socket_to_a_public_address_is_refused():
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        with pytest.raises(NetworkBlockedInTests):
            s.connect(("8.8.8.8", 53))
    finally:
        s.close()


def test_loopback_still_works():
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind(("127.0.0.1", 0))
    server.listen(1)
    port = server.getsockname()[1]
    accepted = []
    t = threading.Thread(target=lambda: accepted.append(server.accept()[0]))
    t.start()
    client = socket.create_connection(("127.0.0.1", port), timeout=5)
    t.join(timeout=5)
    client.close()
    for conn in accepted:
        conn.close()
    server.close()
    assert accepted


def test_the_dataset_listing_does_not_download():
    """The test that used to spend ~19 minutes downloading every cohort:
    listing datasets must finish without any outbound call succeeding."""
    import time

    from fastapi.testclient import TestClient

    from app.api import main as api_main

    api_main._get_dataset.cache_clear()
    start = time.monotonic()
    response = TestClient(api_main.app).get("/datasets")
    assert response.status_code == 200
    assert time.monotonic() - start < 60

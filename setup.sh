mkdir -p ~/.streamlit

cat > ~/.streamlit/config.toml << EOF
[server]
headless = true
enableCORS = false
port = $PORT
EOF

echo "$STREAMLIT_SECRETS_B64" | base64 -d > ~/.streamlit/secrets.toml

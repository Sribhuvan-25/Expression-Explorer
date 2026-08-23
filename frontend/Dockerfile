FROM node:22-slim AS build
WORKDIR /app
COPY package*.json ./
RUN npm install
COPY . .
# API_BASE_URL is baked in at build time since Vite env vars are static —
# point this at wherever the backend actually runs in this deployment.
ARG API_BASE_URL=/api
ENV VITE_API_BASE_URL=$API_BASE_URL
RUN npm run build

FROM nginx:alpine
COPY --from=build /app/dist /usr/share/nginx/html
COPY nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 80

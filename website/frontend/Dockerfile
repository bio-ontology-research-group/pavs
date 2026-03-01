FROM node:20-alpine AS builder

# Leave blank to use relative paths (nginx proxies /api/ to backend).
# Set to the full public URL when the API is on a different origin.
ARG VITE_API_URL=
ENV VITE_API_URL=${VITE_API_URL}

WORKDIR /app
COPY package.json ./
RUN npm install
COPY . .
RUN npm run build

FROM nginx:alpine
COPY --from=builder /app/build /usr/share/nginx/html
COPY nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
